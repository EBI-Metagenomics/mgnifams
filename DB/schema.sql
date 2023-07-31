CREATE TABLE Proteins (
    ID VARCHAR(16) PRIMARY KEY,
    Sequence TEXT NOT NULL
);

CREATE TABLE Families (
    ID VARCHAR(16) PRIMARY KEY,
    IsKnown BOOLEAN,
    MSA BLOB,
    HMM BLOB,
    FOREIGN KEY(ID) REFERENCES Proteins(ID)
);

CREATE TABLE Clustering (
    ProteinID VARCHAR(16),
    FamilyID VARCHAR(16),
    PRIMARY KEY (ProteinID, FamilyID),
    FOREIGN KEY(ProteinID) REFERENCES Proteins(ID)
    FOREIGN KEY(FamilyID) REFERENCES Families(ID)
);

CREATE TABLE Annotations (
    ID INTEGER PRIMARY KEY AUTOINCREMENT,
    FamilyID VARCHAR(16),
    Annotation VARCHAR(30) NOT NULL,
    Description TEXT,
    Source VARCHAR(30),
    IsKnown BOOLEAN,
    FOREIGN KEY(FamilyID) REFERENCES Families(ID)
);