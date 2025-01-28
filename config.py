
class DefaultConfig(object):

    acid_one_hot = [0 for i in range(20)]
    acid_idex = {j:i for i,j in enumerate("ACDEFGHIKLMNPQRSTVWY")} #ACDEFGHIKLMNPQRSTVWYX

    max_sequence_length = 500

parameter = {
    'win_size':13,
}
# [PeatureTool]
Psi_Blast_Path = "./ncbi-blast-2.13.0+/"

# [Database]
DB_PATH ="./nr/ "

# [Result]
Result_Path = ".MMELF-PPI/out/"
