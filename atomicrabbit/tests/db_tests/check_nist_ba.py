"""Check what NIST returns for Ba II to understand term_symbol format."""
from atomicrabbit.db.nist.load_tables import parse_level
from astroquery.nist import Nist
import astropy.units as u

table = Nist.query(400 * u.nm, 700 * u.nm, linename='Ba II')

seen = set()
print("Configuration | term_symbol | J  (from NIST parse_level)")
print("-" * 60)
for row in table:
    lower = parse_level(str(row['Lower level']), 0)
    upper = parse_level(str(row['Upper level']), 0)
    for lvl in [lower, upper]:
        if lvl is None:
            continue
        key = (lvl['configuration'], lvl['term_symbol'], lvl['J'])
        if key not in seen:
            seen.add(key)
            print(f"  {lvl['configuration']!r:20}  {lvl['term_symbol']!r:12}  J={lvl['J']}")
