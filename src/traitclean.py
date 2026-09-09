"""Normalisation of the epidemiological grouping variables.

The GenBank-derived metadata mixes synonyms and granularity levels within a
single column.  Left uncleaned, this fabricates epidemiological differences that
have nothing to do with the inference method:

  host      'bovine' and 'cattle' are the same species, as are 'ovine'/'sheep'
            and 'porcine'/'pig'.  Uncleaned, a cattle-to-cattle transmission is
            scored as a cross-species jump whenever the two records happened to
            be annotated by different submitters.
  division  SARS-CoV-2 USA carries 'WA', 'Washington', 'Washington,King County',
            'Washington, Yakima County' and 'King County,WA' for what is one
            state; H5N1 carries both 'TX'/'Texas' and 'OH'/'Ohio'.  Uncleaned,
            the introduction count is dominated by annotation style.
  division  EBOV Sierra Leone mixes the district 'western area' with its two
            sub-districts 'western rural' and 'western urban', while the other
            three levels are districts.  Mixed granularity makes the count of
            "introductions per district" meaningless, so all three are mapped to
            the district 'Western Area'.

Every mapping below is explicit and auditable; nothing is inferred at runtime.
`normalise()` reports how many labels it changed so the effect is visible in the
results rather than hidden.
"""
import collections

HOST = {
    "bovine": "cattle", "cattle": "cattle", "cow": "cattle", "dairy cow": "cattle",
    "ovine": "sheep", "sheep": "sheep",
    "porcine": "pig", "swine": "pig", "pig": "pig",
    "homo sapiens": "human", "human": "human",
}

US_STATE = {
    "al": "AL", "alabama": "AL", "ak": "AK", "alaska": "AK", "az": "AZ", "arizona": "AZ",
    "ar": "AR", "arkansas": "AR", "ca": "CA", "california": "CA", "co": "CO",
    "colorado": "CO", "ct": "CT", "connecticut": "CT", "de": "DE", "delaware": "DE",
    "dc": "DC", "district of columbia": "DC", "fl": "FL", "florida": "FL",
    "ga": "GA", "georgia": "GA", "hi": "HI", "hawaii": "HI", "id": "ID", "idaho": "ID",
    "il": "IL", "illinois": "IL", "in": "IN", "indiana": "IN", "ia": "IA", "iowa": "IA",
    "ks": "KS", "kansas": "KS", "ky": "KY", "kentucky": "KY", "la": "LA",
    "louisiana": "LA", "me": "ME", "maine": "ME", "md": "MD", "maryland": "MD",
    "ma": "MA", "massachusetts": "MA", "mi": "MI", "michigan": "MI", "mn": "MN",
    "minnesota": "MN", "ms": "MS", "mississippi": "MS", "mo": "MO", "missouri": "MO",
    "mt": "MT", "montana": "MT", "ne": "NE", "nebraska": "NE", "nv": "NV",
    "nevada": "NV", "nh": "NH", "new hampshire": "NH", "nj": "NJ", "new jersey": "NJ",
    "nm": "NM", "new mexico": "NM", "ny": "NY", "new york": "NY", "nc": "NC",
    "north carolina": "NC", "nd": "ND", "north dakota": "ND", "oh": "OH", "ohio": "OH",
    "ok": "OK", "oklahoma": "OK", "or": "OR", "oregon": "OR", "pa": "PA",
    "pennsylvania": "PA", "ri": "RI", "rhode island": "RI", "sc": "SC",
    "south carolina": "SC", "sd": "SD", "south dakota": "SD", "tn": "TN",
    "tennessee": "TN", "tx": "TX", "texas": "TX", "ut": "UT", "utah": "UT",
    "vt": "VT", "vermont": "VT", "va": "VA", "virginia": "VA", "wa": "WA",
    "washington": "WA", "wv": "WV", "west virginia": "WV", "wi": "WI",
    "wisconsin": "WI", "wy": "WY", "wyoming": "WY",
}

SL_DISTRICT = {
    "western area": "Western Area", "western rural": "Western Area",
    "western urban": "Western Area", "port loko": "Port Loko",
    "bombali": "Bombali", "kambia": "Kambia",
}


def _us_state(v):
    """'Washington, Yakima County' / 'King County,WA' / 'WA' -> 'WA'."""
    parts = [p.strip() for p in v.replace(";", ",").split(",") if p.strip()]
    for p in parts:                       # any token that is a state name/code wins
        k = p.lower()
        if k in US_STATE:
            return US_STATE[k]
    return v.strip()


# which cleaner applies to which (dataset, column)
RULES = {
    ("fmdv_uk_2001", "host"): lambda v: HOST.get(v.strip().lower(), v.strip().lower()),
    ("fmdv_uk_2007", "host"): lambda v: HOST.get(v.strip().lower(), v.strip().lower()),
    ("sarscov2_usa_early2020", "division"): _us_state,
    ("h5n1_dairy_cattle_2024", "division"): _us_state,
    ("h5n1_dairy_cattle_2024", "host"): lambda v: HOST.get(v.strip().lower(), v.strip().lower()),
    ("ebov_sierraleone_2014", "division"):
        lambda v: SL_DISTRICT.get(v.strip().lower(), v.strip()),
}


def normalise(dsid, col, trait):
    """Return (cleaned_trait, report).  `trait` is {taxon_label -> raw value}."""
    fn = RULES.get((dsid, col))
    if fn is None:
        return dict(trait), dict(applied=False, n_levels_before=len(set(trait.values())),
                                 n_levels_after=len(set(trait.values())), n_changed=0)
    out = {k: fn(v) for k, v in trait.items()}
    changed = sum(1 for k in trait if out[k] != trait[k])
    return out, dict(applied=True,
                     n_levels_before=len(set(trait.values())),
                     n_levels_after=len(set(out.values())),
                     n_changed=changed,
                     levels_after=dict(sorted(collections.Counter(out.values()).items(),
                                              key=lambda kv: -kv[1])))
