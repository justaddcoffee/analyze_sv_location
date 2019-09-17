from pbga import H2DbManager
import logging

import functools
import gzip

LOG = logging.getLogger(__name__)
LOG.setLevel("INFO")

MIN_OBS = 10  # minimum times to observe SV in order to trust frequency value
MIN_FREQ = 0.05  # minimum frequency for SV to be considered benign

# these are arbitrary scores:
SCORE_BENIGN = 0
SCORE_PATHOGENIC = 1000

config = {
    'DBFILE': "/Users/jtr4v/projects/MI/sv_project/pbga/pbga-main/data/hg19_sv_database.mv.db",
    'USER': "sa",
    'PW': "sa",
    'outfile': 'sv_track.bed.gz',
    'queries':
        {
            'clingen_haploinsuffiency': {
                'pathogenic': f"select CONTIG, START, END, CONCAT('clingenHI gene:', GENE_SYMBOL), {SCORE_PATHOGENIC} "
                              "from PBGA.CLINGEN_HAPLOINSUFFICIENCY",
            },
            'dbvar_variants': {
                'benign': f"SELECT CHR_ONE, POS_ONE, POS_TWO, CONCAT('dbvar ', DBVAR_ACC, ' ', SV_TYPE), {SCORE_BENIGN} "
                          "FROM PBGA.DBVAR_VARIANTS "
                          "WHERE CLNSIG = 'BENIGN' "
                          "OR CLNSIG = 'LIKELY_BENIGN' "
                          f"OR (ALLELE_FREQUENCY != 'NaN' AND ALLELE_FREQUENCY > {MIN_FREQ} AND ALLELE_COUNT > {MIN_OBS})",
                'pathogenic': f"SELECT CHR_ONE, POS_ONE, POS_TWO, CONCAT('dbvar ', DBVAR_ACC, ' ', SV_TYPE), {SCORE_PATHOGENIC}  "
                              "FROM PBGA.DBVAR_VARIANTS "
                              "WHERE CLNSIG = 'PATHOGENIC' "
                              "OR CLNSIG = 'LIKELY_PATHOGENIC' "
            },
            'decipher_cnv': {
                'benign': f"SELECT CONTIG, POS_ONE, POS_TWO, CONCAT('decipher ', POPULATION_CNV_ID, ' ', SV_TYPE), {SCORE_BENIGN} "
                          "FROM PBGA.DECIPHER_CNV "
                          "WHERE FREQUENCY != 'NaN' "
                          f"AND FREQUENCY > {MIN_FREQ} "
                          f"AND OBSERVATIONS > {MIN_OBS}"
            },
            'dgv_variant': {
                'benign': f"SELECT CONTIG, POS_ONE, POS_TWO, CONCAT('dgv ', ACCESSION, ' ', SV_TYPE, ':', DGV_SUBTYPE), {SCORE_BENIGN} "
                          "FROM PBGA.DGV_VARIANTS "
                          f"WHERE ((OBSERVED_GAINS + OBSERVED_LOSSES)/SAMPLE_SIZE > {MIN_FREQ}) "
                          f"AND SAMPLE_SIZE > {MIN_OBS}"
            },
            'gnomad_sv': {
                'benign': f"SELECT CHR_ONE, POS_ONE, POS_TWO, CONCAT('gnomad ', ID, ' ', SV_TYPE), {SCORE_BENIGN} "
                          "FROM PBGA.GNOMAD_SV "
                          f"WHERE POPMAX_AF > {MIN_FREQ} "
                          f"AND AC > {MIN_OBS}"
            },
            'gonl': {
                'benign': f"SELECT CHR_ONE, POS_ONE, POS_TWO, CONCAT('gonl ', ID, '/', DGV_ID, ' ', SV_TYPE), {SCORE_BENIGN} "
                          "FROM PBGA.GONL "
                          f"WHERE AF > {MIN_FREQ} "
                          f"AND AN > {MIN_OBS}"
            },
            'haploinsufficiency': {
                'benign': f"select CONTIG, START, END, CONCAT('HI gene:', GENE_SYMBOL), {SCORE_PATHOGENIC} "
                                  "from PBGA.HAPLOINSUFFICIENCY"
            },
            'dbvar_variants': {
                'benign': f"SELECT CONTIG, POS_ONE, POS_TWO, CONCAT('isca ', ID, ' ', SV_TYPE), {SCORE_BENIGN} "
                          "FROM PBGA.ISCA "
                          "WHERE CLNSIG = 'BENIGN' "
                          "OR CLNSIG = 'LIKELY_BENIGN' ",
                'pathogenic': f"SELECT CONTIG, POS_ONE, POS_TWO, CONCAT('isca ', ID, ' ', SV_TYPE), {SCORE_PATHOGENIC} "
                              "FROM PBGA.ISCA "
                              "WHERE CLNSIG = 'PATHOGENIC' "
                              "OR CLNSIG = 'LIKELY_PATHOGENIC' "
            }
    }
}


def row_to_bed_line(x1, x2):
    return "\t".join([str(x1), str(x2)])


with H2DbManager(db_path=config.get('DBFILE'),
                 user=config.get('USER'),
                 password=config.get('PW')) as h2:
    with h2.get_connection() as conn:

        with gzip.open(config.get("outfile"), 'wt') as f:

            queries = config.get("queries")
            for query in queries:

                for qtype in queries.get(query):
                    LOG.warning("querying table {} for {} SVs".format(query, qtype))
                    select = queries.get(query).get(qtype)

                    with conn.cursor() as cursor:
                        cursor.execute(select)
                        for row in cursor:
                            f.write((functools.reduce(row_to_bed_line, row) + "\n"))
