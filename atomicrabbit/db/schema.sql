-- Spectroscopy database schema for atom-laser interaction simulations
-- SQLite dialect

PRAGMA foreign_keys = ON;

-- ---------------------------------------------------------------------
-- Atoms / isotopes
-- ---------------------------------------------------------------------
CREATE TABLE IF NOT EXISTS atoms (
    atom_id       INTEGER PRIMARY KEY,
    element       TEXT NOT NULL,          -- e.g. 'Rb'
    mass_number   INTEGER,                -- e.g. 87  (NULL if not isotope-specific)
    nuclear_spin  REAL,                   -- I
    mu_I          REAL,                   -- nuclear magnetic dipole moment, in nuclear magnetons
    abundance     REAL,                   -- natural abundance, 0-1, optional
    properties    TEXT,                   -- JSON overflow for novel parameters (e.g. quadrupole moment, mu_I uncertainty)
    UNIQUE (element, mass_number)
);

-- ---------------------------------------------------------------------
-- Papers
-- ---------------------------------------------------------------------
CREATE TABLE IF NOT EXISTS papers (
    paper_id  INTEGER PRIMARY KEY,
    title     TEXT NOT NULL,
    authors   TEXT,
    year      INTEGER,
    doi       TEXT UNIQUE,
    journal   TEXT,
    notes     TEXT
);

-- what kind(s) of data a paper contributes
CREATE TABLE IF NOT EXISTS paper_data_types (
    paper_id   INTEGER NOT NULL REFERENCES papers(paper_id) ON DELETE CASCADE,
    data_type  TEXT NOT NULL CHECK (data_type IN (
                   'energy_level', 'transition', 'isotope_shift',
                   'hyperfine_constant', 'branching_ratio', 'lifetime', 'other'
               )),
    PRIMARY KEY (paper_id, data_type)
);

-- citation graph: citing_paper cites cited_paper
CREATE TABLE IF NOT EXISTS paper_references (
    citing_paper_id  INTEGER NOT NULL REFERENCES papers(paper_id) ON DELETE CASCADE,
    cited_paper_id   INTEGER NOT NULL REFERENCES papers(paper_id) ON DELETE CASCADE,
    PRIMARY KEY (citing_paper_id, cited_paper_id),
    CHECK (citing_paper_id != cited_paper_id)
);

-- ---------------------------------------------------------------------
-- Energy levels
-- ---------------------------------------------------------------------
CREATE TABLE IF NOT EXISTS energy_levels (
    level_id            INTEGER PRIMARY KEY,
    atom_id             INTEGER NOT NULL REFERENCES atoms(atom_id),
    configuration       TEXT,             -- e.g. '5s 5p'
    term_symbol         TEXT,             -- e.g. '2P3/2'
    L                   REAL,             -- orbital angular momentum quantum number
    S                   REAL,             -- spin angular momentum quantum number
    J                   REAL,
    F                   REAL,             -- hyperfine level, nullable if fine-structure only
    parity              TEXT,
    energy_cm1          REAL NOT NULL,    -- relative to ground state, cm^-1
    energy_uncertainty  REAL,
    hyperfine_A_MHz     REAL,
    hyperfine_B_MHz     REAL,
    lifetime_s          REAL,
    is_recommended      INTEGER NOT NULL DEFAULT 0,  -- 1 = preferred value for simulation queries
    properties          TEXT,                         -- JSON overflow for novel parameters
    source_type         TEXT NOT NULL CHECK (source_type IN ('nist', 'paper')),
    source_paper_id     INTEGER REFERENCES papers(paper_id),
    CHECK (
        (source_type = 'paper' AND source_paper_id IS NOT NULL) OR
        (source_type = 'nist'  AND source_paper_id IS NULL)
    )
);

CREATE INDEX IF NOT EXISTS idx_levels_atom ON energy_levels(atom_id);

-- ---------------------------------------------------------------------
-- Transitions
-- ---------------------------------------------------------------------
CREATE TABLE IF NOT EXISTS transitions (
    transition_id          INTEGER PRIMARY KEY,
    atom_id                INTEGER NOT NULL REFERENCES atoms(atom_id),
    upper_level_id         INTEGER NOT NULL REFERENCES energy_levels(level_id),
    lower_level_id         INTEGER NOT NULL REFERENCES energy_levels(level_id),
    wavelength_nm          REAL,
    frequency_THz          REAL,
    einstein_A_s1          REAL,          -- decay rate Gamma, s^-1
    linewidth_MHz          REAL,          -- natural linewidth (= A / 2pi for single-channel decay)
    isotope_shift_MHz      REAL,
    isotope_shift_ref_mass INTEGER,       -- isotope the shift is quoted relative to
    branching_ratio        REAL,
    is_recommended         INTEGER NOT NULL DEFAULT 0,  -- 1 = preferred value for simulation queries
    properties             TEXT,                         -- JSON overflow for novel parameters
    source_type            TEXT NOT NULL CHECK (source_type IN ('nist', 'paper')),
    source_paper_id        INTEGER REFERENCES papers(paper_id),
    CHECK (
        (source_type = 'paper' AND source_paper_id IS NOT NULL) OR
        (source_type = 'nist'  AND source_paper_id IS NULL)
    ),
    UNIQUE (upper_level_id, lower_level_id, source_type, source_paper_id)
);

CREATE INDEX IF NOT EXISTS idx_transitions_atom  ON transitions(atom_id);
CREATE INDEX IF NOT EXISTS idx_transitions_upper ON transitions(upper_level_id);
CREATE INDEX IF NOT EXISTS idx_transitions_lower ON transitions(lower_level_id);
