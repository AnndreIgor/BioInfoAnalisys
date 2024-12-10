import random
import os

# match algoritmo:
#     case 'clustalo':
#         params = {
#             'outfmt': 'clu',
#         }
#         tags = ['auto']
#     case 'mafft':
#         params = {
#             'retree':2
#         }
#         tags = ['auto']
#     case 't_coffee':
#         params = {
#             'output':'clustalw'
#         }
#     case 'muscle':
#         tags = ['clw']

def __sortear(params: dict, tags: list):
    # Sorteando os parametros do dicionário
    num_chaves = random.randint(0, len(params.keys()))
    chaves_aleatorias = random.sample(list(params.keys()), num_chaves)
    _params = {chave: params[chave] for chave in chaves_aleatorias}
    _tags = random.sample(tags,  random.randint(0, len(tags)))

    return _params, _tags


def sort_params(algoritmo: str) -> tuple:
    match algoritmo:
        case 'clustalo':
            params = {
                'cluster-size': random.randint(1, 10),
                'trans': random.randint(1, 3),
                'threads': random.randint(1, os.cpu_count())
            }

            tags = ('dealign', 'is-profile', 'pileup', 'full', 'full-iter', 'use-kimura', 'percent-id', 'residuenumber', 'auto')
            params, tags = __sortear(params, tags)
            params['outfmt'] = 'clu' # Esse paramtro é obrigatório
        
        case 'mafft': 
            params = {
                'op': random.randint(1, 3),
                'ep': random.randint(1, 3),
                'maxiterate': random.randint(10, 20),
            }

            tags = ('reorder', 'quiet', 'dash')
            params, tags = __sortear(params, tags)

            tags.append('clustalout')
        
        case 'muscle':
            params = {
                'maxiters': random.randint(2, 20),
                'maxhours': random.randint(1, 2)
            }
            tags = ('quiet', 'diags')
            params, tags = __sortear(params, tags)
            tags.append('clw')

        case 'clustalw':
            params = {
                'SEQNOS': random.choice(['OFF', 'ON']),
            }

            tags = ('QUICKTREE', 'NEGATIVE', 'QUIET')
            params, tags = __sortear(params, tags)
            params['OUTPUT'] = 'CLUSTAL' # Esse parametro é obrigatório

        case 'probcons':
            params = {
                'c': random.randint(0, 5),
                'ir': random.randint(0, 1000),
                'pre': random.randint(0, 20),

            }
            tags = ('pairs', 'viterbi', 'v')

            params, tags = __sortear(params, tags)
            tags.append('clustalw')

        case 't_coffee':
            params = {

            }
            tags = ()

            params, tags = __sortear(params, tags)
            params['output'] = 'clustalw' # Esse parametro é obrigatório

    return params, tags

if __name__ == '__main__':
    params, tags = sort_params('mafft')
    print(params)
    print(tags)