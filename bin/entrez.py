#!/usr/bin/env python3
#
# A polite Entrez client, shared by resolve_reference.py and build_reference.py.
# Every call carries tool= and email= as NCBI's usage policy requires, requests
# are spaced to the rate limit (10/s with an API key in NCBI_API_KEY, 3/s
# without), and a throttle or an error is raised rather than returned empty:
# an empty hit set must mean "this taxon has no sequences", never "NCBI said no".
# Responses are cached on disk so a re-run costs nothing.
#
import hashlib, json, os, sys, time, urllib.error, urllib.parse, urllib.request

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
TOOL = "nectarflow-lowe2022"


class Throttled(Exception):
    pass


class Entrez:
    #
    # email is required by the usage policy; it is read from --email or from
    # NCBI_EMAIL. cache is a directory of response bodies keyed by the request,
    # with the API key left out of the key so a cache survives getting one
    def __init__(self, email=None, api_key=None, cache=None, tool=TOOL, sleep=None):
        self.email = email or os.environ.get("NCBI_EMAIL", "")
        self.api_key = api_key or os.environ.get("NCBI_API_KEY", "")
        if not self.email:
            raise SystemExit("entrez: an email address is required (--email or NCBI_EMAIL)")
        self.tool, self.cache = tool, cache
        self.sleep = sleep if sleep is not None else (0.11 if self.api_key else 0.35)
        self.last = 0.0
        self.n_req = self.n_hit = 0
        if cache:
            os.makedirs(cache, exist_ok=True)

    #
    # Builds the query string, adding the policy fields, and returns it with the
    # cache key computed before the API key is added
    def _url(self, cgi, params):
        p = dict(params, tool=self.tool, email=self.email)
        key = hashlib.sha256((cgi + "?" + urllib.parse.urlencode(sorted(p.items()))).encode()).hexdigest()
        if self.api_key:
            p["api_key"] = self.api_key
        return EUTILS + cgi + "?" + urllib.parse.urlencode(p), key

    #
    # One request, from the cache if it is there. HTTP 429 and the "API rate
    # limit exceeded" body are retried with a growing wait and then raised;
    # everything else is retried a few times and then raised as well
    def get(self, cgi, params, post=False, tries=6, fresh=False):
        url, key = self._url(cgi, params)
        path = os.path.join(self.cache, key) if self.cache else None
        if path and os.path.exists(path) and not fresh:
            self.n_hit += 1
            return open(path, encoding="utf-8").read()
        body = None
        if post:
            url, qs = url.split("?", 1)
            body = qs.encode()
        for i in range(tries):
            wait = self.sleep - (time.time() - self.last)
            if wait > 0:
                time.sleep(wait)
            self.last = time.time()
            self.n_req += 1
            try:
                req = urllib.request.Request(url, data=body)
                with urllib.request.urlopen(req, timeout=300) as r:
                    text = r.read().decode("utf-8", "replace")
                if "API rate limit exceeded" in text[:400]:
                    raise Throttled("rate limit in body")
                if path:
                    # via a temp file in the same directory, then replace: the cache is trusted on
                    # the strength of the file existing, so a fetch interrupted midway through the
                    # write would otherwise leave a truncated body that every later run believes
                    tmp = f"{path}.part{os.getpid()}"
                    with open(tmp, "w", encoding="utf-8") as fh:
                        fh.write(text)
                    os.replace(tmp, path)
                return text
            except urllib.error.HTTPError as e:
                if e.code == 429:
                    time.sleep(2 * (i + 1))
                    continue
                if e.code in (500, 502, 503, 504):
                    time.sleep(2 * (i + 1))
                    continue
                raise
            except (Throttled, urllib.error.URLError, TimeoutError, OSError) as e:
                print(f"entrez retry {i + 1}: {e}", file=sys.stderr)
                time.sleep(2 * (i + 1))
        raise SystemExit(f"entrez: gave up after {tries} tries: {cgi} {params.get('term', params.get('id', ''))[:80]}")

    #
    # esearch with usehistory, returning the hit count and the history handle.
    # The sort order NCBI returns is never relied on: callers page the whole set
    # and sort it themselves, because an unknown sort= is answered HTTP 200 with
    # the default order and no warning
    def esearch(self, db, term, extra=None):
        p = {"db": db, "term": term, "usehistory": "y", "retmax": "0", "retmode": "json"}
        p.update(extra or {})
        d = json.loads(self.get("esearch.fcgi", p), strict=False)["esearchresult"]
        if "ERROR" in d:
            raise SystemExit(f"entrez esearch error: {d['ERROR']} for {term}")
        return int(d["count"]), d.get("webenv", ""), d.get("querykey", "")

    #
    # Pages the whole hit set of a previous esearch as document summaries,
    # yielding one dict per record (accession.version, length, organism, taxid,
    # title). 500 at a time, which is what esummary is comfortable with
    def esummary_all(self, db, webenv, key, count, step=500):
        for start in range(0, count, step):
            p = {"db": db, "query_key": key, "WebEnv": webenv, "retstart": str(start),
                 "retmax": str(min(step, count - start)), "retmode": "json", "version": "2.0"}
            # An empty page is usually NCBI being momentarily unwilling rather than the
            # records being absent, so ask again a few times before giving up on the run.
            for attempt in range(5):
                # strict=False: NCBI leaves raw control characters in some record titles
                d = json.loads(self.get("esummary.fcgi", p, fresh=attempt > 0),
                               strict=False).get("result", {})
                if d.get("uids"):
                    break
                time.sleep(2 ** attempt)
            else:
                raise SystemExit(f"entrez esummary returned nothing at retstart={start} of {count}"
                                 " after 5 tries")
            for uid in d["uids"]:
                r = d[uid]
                yield {"accession": r.get("accessionversion", ""), "length": int(r.get("slen", 0) or 0),
                       "organism": r.get("organism", ""), "taxid": int(r.get("taxid", 0) or 0),
                       "title": r.get("title", "")}

    #
    # Posts a list of accessions and fetches them back as FASTA. epost silently
    # truncates long id lists - 2000 posted has come back as 1756 with no error
    # and an empty InvalidIdList - so the chunk is small and the caller checks
    # that every accession it asked for came back
    def fetch_fasta(self, accessions, chunk=200):
        out = {}
        for i in range(0, len(accessions), chunk):
            part = accessions[i:i + chunk]
            d = self.get("epost.fcgi", {"db": "nuccore", "id": ",".join(part)}, post=True)
            webenv = _between(d, "<WebEnv>", "</WebEnv>")
            key = _between(d, "<QueryKey>", "</QueryKey>")
            if not webenv:
                raise SystemExit("entrez epost returned no WebEnv")
            text = self.get("efetch.fcgi", {"db": "nuccore", "query_key": key, "WebEnv": webenv,
                                            "rettype": "fasta", "retmode": "text"})
            out.update(parse_fasta(text))
        return out


#
# The text between two tags, or an empty string
def _between(s, a, b):
    i = s.find(a)
    return "" if i < 0 else s[i + len(a):s.find(b, i)]


#
# Reads FASTA text into {accession.version: uppercase ungapped sequence}. The
# accession is the first token of the header, which is what efetch writes
def parse_fasta(text):
    out, name, seq = {}, None, []
    for line in text.splitlines():
        if line.startswith(">"):
            if name:
                out[name] = "".join(seq)
            name, seq = line[1:].split()[0], []
        elif line.strip():
            seq.append("".join(c for c in line.strip().upper() if c.isalpha()))
    if name:
        out[name] = "".join(seq)
    return out
