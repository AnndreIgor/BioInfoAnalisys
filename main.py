from Bio import SeqIO
from Bio.Align.Applications import ClustalwCommandline, MafftCommandline, ClustalOmegaCommandline, MuscleCommandline, TCoffeeCommandline, ProbconsCommandline
from Bio import AlignIO, Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.Phylo.Applications import RaxmlCommandline

import os


def sequence_align(input_file: str,
                   output_file: str,
                   algorithm: str,
                   output_dir: str):

    output_file = os.path.join(output_dir, output_file)

    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    lista_algoritmos = [
        'muscle',
        'clustalw',
        'clustalo',
        'mafft',
        'probcons',
        't-coffee'
    ]

    if algorithm not in lista_algoritmos:
        print('Algoritmo de alinhamento ainda não implementado')
        return None

    # Carregar as sequências a serem alinhadas
    sequences = SeqIO.parse(input_file, "fasta")

    # # Escrever as sequências em um arquivo temporário no formato FASTA
    # SeqIO.write(sequences, "temp.fasta", "fasta")

    match algorithm.lower():

        case "muscle":
            # Definir a linha de comando para executar o MUSCLE
            muscle_cline = MuscleCommandline(input=input_file, out=output_file)

            # Executar o MUSCLE
            stdout, stderr = muscle_cline()

        case "clustalw":
            # Definir a linha de comando para executar o ClustalW
            clustalw_cline = ClustalwCommandline("clustalw", infile=input_file, outfile=output_file)

            # Executar o ClustalW
            stdout, stderr = clustalw_cline()

        case "clustalo":
            # Definir a linha de comando para executar o Clustal Omega
            clustalo_cline = ClustalOmegaCommandline(
                infile=input_file,
                outfile=output_file,
                verbose=True,
                auto=True
            )

            # Executar o Clustal Omega
            stdout, stderr = clustalo_cline()

        case "mafft":
            # Definir a linha de comando para executar o MAFFT
            mafft_cline = MafftCommandline(input=input_file)

            # Executar o MAFFT
            stdout, stderr = mafft_cline()
            with open(output_file, "w") as handle:
                handle.write(stdout)

        case "t-coffee":
            # Crie o objeto de linha de comando para o T-Coffee
            tcoffee_cline = TCoffeeCommandline(infile=input_file,
                                               output="clustalw",
                                               outfile=output_file)
            tcoffee_cline()

        case "probcons":
            probcons_cline = ProbconsCommandline(input=input_file, clustalw=True)
            stdout, stderr = probcons_cline()

            with open(output_file, "w") as handle:
                handle.write(stdout)


def construir_arvore_filogenetica(sequencia_alinhada: str,
                                  caminho_saida: str,
                                  formato_saida: str,
                                  metodo: str
                                  ) -> None:
    if formato_saida not in [
        'newick'
    ]:
        print(f'Formato de saída {formato_saida} não aceito')
        return None

    # Carregar o arquivo de alinhamento de sequências
    aln = AlignIO.read(sequencia_alinhada, "clustal")

    match metodo.lower():
        case 'distancetreeconstructor':
            # Calcular as distâncias entre as sequências
            calculator = DistanceCalculator('identity')
            dm = calculator.get_distance(aln)

            # Construir a árvore filogenética
            constructor = DistanceTreeConstructor()
            tree = constructor.upgma(dm)

            # Salvar a árvore em um arquivo
            Phylo.write(tree, os.path.join(caminho_saida, "arvore_filogenetica.nwk"), formato_saida)

        case 'raxml':
            ## Executar o RAxML
            AlignIO.write(aln, "temp.phy", "phylip")

            raxml_cline = RaxmlCommandline(
                sequence="temp.phy",
                model="GTRGAMMA",
                rapid_bootstrap=100,
                num_replicates=10,
                parsimony_seed=12345,
                name="output"
            )
            stdout, stderr = raxml_cline()

            stdout, stderr = raxml_cline()

    print('Arvore construida com sucesso')


def mostrar_arvore(caminho_arvore: str,
                   formato: str):
    # Ler a árvore filogenética de saída
    tree = Phylo.read(caminho_arvore, formato)

    # Visualizar a árvore
    Phylo.draw(tree)

# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    input_file = os.path.join('files', 'input', 'ls_orchid.fasta')
    output_dir = os.path.join('files', 'output', 'sequencias_alinhadas')
    sequence_align(input_file,
                   'aligned.aln',
                   'probcons',
                   output_dir)

    input_folder = os.path.join('files', 'output', 'sequencias_alinhadas')
    output_folder = os.path.join('files', 'output', 'arvores_filogeneticas')
    construir_arvore_filogenetica(
        os.path.join(input_folder, 'aligned.aln'),
        output_folder,
        'newick',
        'distancetreeconstructor'
    )

    caminho_arvore = os.path.join(output_folder, "arvore_filogenetica.nwk")
    mostrar_arvore(caminho_arvore, 'newick')