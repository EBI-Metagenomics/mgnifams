import psycopg2

# PostgreSQL database parameters
database_name = 'mgnproteindev'
user = 'mgnprotein_rw'
password = 'changethis' # TODO read from file outside git after changed
host = 'pgsql-hlvm-085.ebi.ac.uk'

conn = psycopg2.connect(dbname=database_name, user=user, password=password, host=host)
cursor = conn.cursor()
cursor.execute("SELECT mgyp, sequence FROM sequence_explorer_protein")

with open('pgsql.fasta', 'w') as f:
    for row in cursor:
        f.write(f">{row[0]}\n{row[1]}\n")

cursor.close()
conn.close()
