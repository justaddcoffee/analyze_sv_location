from pbga import H2DbManager
import logging
import functools
import gzip
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('-q', '--query_type', type=str,
                    choices=list(["benign", "pathogenic"]), nargs='+', required=True,
                    help='query for benign SVs, pathogenic SVs, or both')
parser.add_argument('-o', '--outfile', type=str, required=True,
                    help='name for output file')
parser.add_argument('-c', '--correct_start_end', action='store_true',
                    help='should we swap start and end if start > end? [false]')
parser.add_argument('-d', '--db', type=str, required=True,
                    help='h2 database file (include .sv.db extension)')
parser.add_argument('-u', '--user', type=str, default='sa', help='db username')
parser.add_argument('-p', '--pw', type=str, default='sa', help='db password')
parser.add_argument('-b', '--min_obs', type=int, default=10, help=
                    'minimum times to observe SV in order to trust frequency value [10]')
parser.add_argument('-f', '--min_freq', type=float, default=0.05, help=
                    'minimum frequency for SV to be considered benign')
parser.add_argument('--score_path', type=int, default=1000, help=
                    'score to output in bed for pathogenic variants')
parser.add_argument('--score_benign', type=int, default=0, help=
                    'score to output in bed for pathogenic variants')

opts = parser.parse_args()
print('query type was:', opts.query_type)

LOG = logging.getLogger(__name__)
LOG.setLevel("INFO")

queries = {
    'clingen_haploinsuffiency': {
        'pathogenic': f"select CONTIG, START, END, CONCAT('clingenHI gene:', GENE_SYMBOL), {opts.score_path} "
                      "from PBGA.CLINGEN_HAPLOINSUFFICIENCY",
    },
    'dbvar_variants': {
        'benign': f"SELECT CHR_ONE, POS_ONE, POS_TWO, CONCAT('dbvar ', DBVAR_ACC, ' ', SV_TYPE), {opts.score_benign} "
                  "FROM PBGA.DBVAR_VARIANTS "
                  "WHERE CHR_ONE = CHR_TWO "
                  "AND (CLNSIG = 'BENIGN' "
                  "OR CLNSIG = 'LIKELY_BENIGN' "
                  f"OR (ALLELE_FREQUENCY != 'NaN' AND ALLELE_FREQUENCY > {opts.min_freq} AND ALLELE_COUNT > {opts.min_obs}))",
        'pathogenic': f"SELECT CHR_ONE, POS_ONE, POS_TWO, CONCAT('dbvar ', DBVAR_ACC, ' ', SV_TYPE), {opts.score_path}  "
                      "FROM PBGA.DBVAR_VARIANTS "
                      "WHERE CLNSIG = 'PATHOGENIC' "
                      "OR CLNSIG = 'LIKELY_PATHOGENIC' "
    },
    'decipher_cnv': {
        'benign': f"SELECT CONTIG, POS_ONE, POS_TWO, CONCAT('decipher ', POPULATION_CNV_ID, ' ', SV_TYPE), {opts.score_benign} "
                  "FROM PBGA.DECIPHER_CNV "
                  "WHERE FREQUENCY != 'NaN' "
                  f"AND FREQUENCY > {opts.min_freq} "
                  f"AND OBSERVATIONS > {opts.min_obs}"
    },
    'dgv_variant': {
        'benign': f"SELECT CONTIG, POS_ONE, POS_TWO, CONCAT('dgv ', ACCESSION, ' ', SV_TYPE, ':', DGV_SUBTYPE), {opts.score_benign} "
                  "FROM PBGA.DGV_VARIANTS "
                  f"WHERE ((OBSERVED_GAINS + OBSERVED_LOSSES)/SAMPLE_SIZE > {opts.min_freq}) "
                  f"AND SAMPLE_SIZE > {opts.min_obs}"
    },
    'gnomad_sv': {
        'benign': f"SELECT CHR_ONE, POS_ONE, POS_TWO, CONCAT('gnomad ', ID, ' ', SV_TYPE), {opts.score_benign} "
                  "FROM PBGA.GNOMAD_SV "
                  "WHERE CHR_ONE = CHR_TWO "
                  f"AND POPMAX_AF > {opts.min_freq} "
                  f"AND AC > {opts.min_obs}"
    },
    'gonl': {
        'benign': f"SELECT CHR_ONE, POS_ONE, POS_TWO, CONCAT('gonl ', ID, '/', DGV_ID, ' ', SV_TYPE), {opts.score_benign} "
                  "FROM PBGA.GONL "
                  "WHERE CHR_ONE = CHR_TWO "
                  f"AND AF > {opts.min_freq} "
                  f"AND AN > {opts.min_obs}"
    },
    'haploinsufficiency': {
        'benign': f"select CONTIG, START, END, CONCAT('HI gene:', GENE_SYMBOL), {opts.score_path} "
                  "from PBGA.HAPLOINSUFFICIENCY"
    },
    'dbvar_variants': {
        'benign': f"SELECT CONTIG, POS_ONE, POS_TWO, CONCAT('isca ', ID, ' ', SV_TYPE), {opts.score_benign} "
                  "FROM PBGA.ISCA "
                  "WHERE CLNSIG = 'BENIGN' "
                  "OR CLNSIG = 'LIKELY_BENIGN' ",
        'pathogenic': f"SELECT CONTIG, POS_ONE, POS_TWO, CONCAT('isca ', ID, ' ', SV_TYPE), {opts.score_path} "
                      "FROM PBGA.ISCA "
                      "WHERE CLNSIG = 'PATHOGENIC' "
                      "OR CLNSIG = 'LIKELY_PATHOGENIC' "
    }
}


def row_to_bed_line(x1, x2):
    return "\t".join([str(x1), str(x2)])


def correct_start_end(row):
    tmp = row[1]
    row[1] = row[2]
    row[2] = tmp
    return row


with H2DbManager(db_path=opts.db, user=opts.user, password=opts.pw) as h2:
    with h2.get_connection() as conn:

        with gzip.open(opts.outfile, 'wt') as f:

            for query in queries:

                for qtype in opts.query_type:
                    LOG.info("querying table {} for {} SVs".format(query, qtype))
                    select = queries.get(query).get(qtype)
                    if queries.get(query).get(qtype) is None:
                        continue

                    with conn.cursor() as cursor:
                        cursor.execute(select)
                        for row in cursor:
                            if row[1] > row[2]:
                                LOG.warning("start is after stop in record: {}".
                                            format(row))
                                if opts.correct_start_end:
                                    LOG.warning("correcting start and end")
                                    row = correct_start_end(list(row))
                            f.write((functools.reduce(row_to_bed_line, row) + "\n"))
