# Attribution and licensing scope

Project software and website: Copyright 2026 Zijia Wang, Haoliang Hong, Xuhui Li, Jingyang Bai, and Binjian Nie.
Manuscript authors: Zijia Wang, Haoliang Hong, Xuhui Li, Jingyang Bai, and Binjian Nie.
This release does not declare a blanket license for the project-authored source, paper or results. Making the
repository public does not change the upstream terms below.

## Isotherm data in the explorer

The interactive explorer (`docs/explorer/data/`) contains measured points of 411 pure-gas isotherm records obtained
from the public API of the NIST/ARPA-E Database of Novel and Emerging Adsorbent Materials (NIST Standard Reference
Database 205, https://adsorption.nist.gov/isodb/, retrieved 2026-11-09), together with quantities that this project
computed from them (compiled envelope tables, certified bands, least-squares fits, lookahead advice). Each record
carries its ISODB filename and the DOI of the article it was digitized from; `scripts/demo_sources.json` lists the
selection. The records are included as small factual excerpts for illustration and verification. They remain subject
to the terms of the database and of the original publications; no license over them is granted or implied here, and
no NIST endorsement is implied. A record will be removed on request of its source.

The Python package does not bundle isotherm data. `isocert.isodb.fetch_isotherm` downloads records from the public
API on the user's request; acquire and use such data under the original terms and attribution requirements.

## Other material

`scripts/render_demos.py` and the page layout follow the authors' companion project
(https://github.com/OliverDOU776/Transporting-Root-Evidence-for-Repeated-Adsorption-Ranking). The project page is an
original implementation; it borrows organizational ideas from academic project pages rather than copying their source
or media.
