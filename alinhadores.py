from Bio import AlignIO

def make_clustalw(file: str, output: str, *args, **kwargs) -> list:
    """Monta o comando para execução do clustalw em linha de comando
    gera 2 arquivos de saida por padrão

    Nesse caso:
    - ORTHOMCL256.aln: Contendo a sequência de ORTHOMCL256 alinhada em formato clustal
    - ORTHOMCL256.dnd: Contendo informações sobre o agrupamento hierárquico das sequências alinhadas.

    Mais informações sobre aplicações biopython e clustalwcommandline: https://biopython.org/docs/1.76/api/Bio.Align.Applications.html

    Args:
        file (str): Arquivo de sequencias .fasta
        output (str): Arquivo de saída .aln

    Returns:
        list: lista de comandos com os parâmetros utilizados
    """

    comand = ['clustalw', f'-INFILE={file}', f'-OUTFILE={output}']

    for key, value in kwargs.items():
        comand.append('-' + key + '=' + str(value))

    for x in args:
        comand.append('-' + x)

    return comand


def make_muscle(file: str, output: str, *args, **kwargs):
    """Cria a lista para executar o muscle via subprocess
    https://drive5.com/muscle5/manual/cmd_align.html

    Args:
        file (str): Arquivo de sequencias .fasta
        output (str): Arquivo de saída .aln

    Returns:
        _type_: lista de comandos
    """
    
    comand = ['muscle', '-in', file, '-out', output]

    for key, value in kwargs.items():
        comand.append('-' + key + ' ' + str(value))

    for x in args:
        comand.append('-' + x)

    return comand


def make_clustalo(file: str, output: str, *args, **kwargs):
    """ http://www.clustal.org/omega/README

    Args:
        file (str): _description_
        output (str): _description_

    Returns:
        _type_: _description_
    """
    comand = ['clustalo', '-i',  file, '-o', output]

    for key, value in kwargs.items():
        comand.append('--' + key + '=' + str(value))

    for x in args:
        comand.append('--' + x)

    return comand


def make_mafft(file: str, *args, **kwargs):
    """_summary_
    https://mafft.cbrc.jp/alignment/software/manual/manual.html

    Args:
        file (str): _description_
        output (str): _description_

    Returns:
        _type_: _description_
    """

    comand = ['mafft']

    for key, value in kwargs.items():
        comand.append('--' + key)
        comand.append(str(value))

    for x in args:
        comand.append('--' + x)

    comand.append(file)

    return comand

def make_probcons(file: str, *args, **kwargs):
    """
    http://probcons.stanford.edu/manual.pdf
    
    """
    comand = ['probcons']

    for key, value in kwargs.items():
        comand.append('-' + key)
        comand.append(str(value))

    for x in args:
        comand.append('-' + x)

    comand.append(file)

    return comand


def make_t_coffee(file: str, outfile: str, *args, **kwargs):
    comand = ['t_coffee', file, '-outfile', outfile]

    for key, value in kwargs.items():
        comand.append('-' + key)
        comand.append(str(value))

    for x in args:
        comand.append('-' + x)

    return comand


def to_clustalw(file: str):
    # Lê o alinhamento em formato FASTA
    alignment = AlignIO.read(file, "fasta")

    # Escreve o alinhamento no formato Clustal
    AlignIO.write(alignment, file, "clustal")

if __name__ == '__main__':
    entrada = "data/full_dataset_plasmodium/PLASMODIUM0.fasta"
    
    saida = "data/out/tmp/PLASMODIUM0_muscle_.aln"
    comand = make_muscle(entrada, saida, 'clw')
    print(" ".join(comand))
    print()

    saida = "data/out/tmp/PLASMODIUM0_clustalw_.aln"
    comand = make_clustalw(entrada, saida, 'QUIET', 'NOPGAP', STATS='stats_file.txt')
    print(" ".join(comand))
    print()

    saida = "data/out/tmp/PLASMODIUM0_clustalo_.aln"
    comand = make_clustalo(entrada, saida, 'auto', threads=4, log='log_file.txt', outfmt='clu')
    print(" ".join(comand))
    print()

    saida = "data/out/tmp/PLASMODIUM0_mafft_.aln"
    comand = make_mafft(entrada, retree=2)
    print(" ".join(comand))
    print()

    saida = "data/out/tmp/PLASMODIUM0_probcons_.aln"
    comand = make_probcons(entrada)
    print(" ".join(comand))
    print()

    saida = "data/out/tmp/PLASMODIUM0_t_coffee_.aln"
    comand = make_t_coffee(entrada, saida, output='clustalw')
    print(" ".join(comand))
    print()


    