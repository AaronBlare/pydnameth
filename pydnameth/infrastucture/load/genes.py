from pydnameth.infrastucture.path import get_data_base_path
import numpy as np
import os.path
import pickle
from tqdm import tqdm
from pydnameth.infrastucture.load.betas import get_line_list


def load_genes(config):
    fn_dict = get_data_base_path(config) + '/' + 'gene_dict.pkl'

    suffix = ''
    if bool(config.experiment.data_params):
        suffix += '_' + str(config.experiment.get_data_params_str())

    fn_txt = get_data_base_path(config) + '/' + 'betas' + suffix + '.txt'
    fn_npz = get_data_base_path(config) + '/' + 'gene' + suffix + '.npz'

    if os.path.isfile(fn_dict) and os.path.isfile(fn_npz):

        f = open(fn_dict, 'rb')
        config.gene_dict = pickle.load(f)
        f.close()

        data = np.load(fn_npz)
        config.gene_data = data['data']

    else:

        config.gene_dict = {}
        config.gene_list = []

        for gene_id, gene in enumerate(config.gene_cpg_dict):
            config.gene_dict[gene] = gene_id
            config.gene_list.append(gene)

        f = open(fn_dict, 'wb')
        pickle.dump(config.gene_dict, f, pickle.HIGHEST_PROTOCOL)
        f.close()

        f = open(fn_txt)
        header_line = f.readline()
        headers = header_line.split('\t')
        headers = [x.rstrip() for x in headers]
        subjects = headers[1:len(headers)]

        config.gene_data = np.zeros((len(config.gene_list), len(subjects)), dtype=np.float32)

        gene_id = 0
        for line in tqdm(f, mininterval=60.0, desc='gene_data creating'):
            line_list = get_line_list(line)
            curr_cpg_name = line_list[0]
            curr_data = list(map(np.float32, line_list[1::]))

            for gene in config.cpg_gene_dict[curr_cpg_name]:
                line_num = config.gene_dict[gene]
                config.gene_data[line_num] += curr_data

            gene_id += 1
        f.close()

        for row_id, row in enumerate(config.gene_data):
            config.gene_data[row_id] /= len(config.gene_cpg_dict[config.gene_list[row_id]])

        np.savez_compressed(fn_npz, data=config.gene_data)
