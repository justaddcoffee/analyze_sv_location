# for the following to work, you need to roll back JPype to 0.6.3:
# `pip install JPype1==0.6.3 --force-reinstall`
# per this SO exchange:
# https://github.com/baztian/jaydebeapi/issues/99

import jaydebeapi
h2_db_prefix = 'jdbc:h2:file:'
h2_location = '/Users/jtr4v/projects/MI/sv_project/pbga/pbga-main/data/hg19/data/sv_database'
h2_jar_file = '/Users/jtr4v/projects/MI/sv_project/exomiser/exomiser-main/exomiser-cli/target/lib/h2-1.4.197.jar'
user = 'sa'
pw = ''
conn = jaydebeapi.connect(jclassname='org.h2.Driver',
                          url=h2_db_prefix + h2_location,
                          driver_args=[user, pw],
                          jars=h2_jar_file)
curs = conn.cursor()
curs.execute("SELECT * FROM PBGA.DBVAR LIMIT 10")
rows = curs.fetchall()
for row in rows:
    print(row)
