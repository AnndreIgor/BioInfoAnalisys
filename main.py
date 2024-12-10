# %%
import os
import subprocess
from pathlib import Path
import pandas as pd
import psutil
import time
from sqlalchemy import create_engine
import random
import logging
from collections import Counter

# %%
from tabelas import *
from metricas import *
from alinhadores import *
from parametros_algoritmos import sort_params

# %%
from Bio import AlignIO, Phylo, SeqIO
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from dendropy import Tree 


logging.basicConfig(
    level=logging.DEBUG,  # Define o nível de log
    format="%(asctime)s [%(levelname)s] %(message)s",  # Formato da mensagem
    datefmt="%Y-%m-%d %H:%M:%S",  # Formato da data
    handlers=[
        logging.FileHandler("app.log"),  # Log em arquivo
    ]
)

# %%
def duplicate_names(file_path: str) -> bool:
    """Verifica se existem sequências duplicadas em um arquivo

    Args:
        file_path (str): Caminho do arquivo

    Returns:
        bool: True se houver sequências duplicadas e False se não houver
    """
    
    names_set = set()

    try:
        for record in SeqIO.parse(file_path, 'fasta'):
            if record.id in names_set:
                return True
            else:
                names_set.add(record.id)

    except FileNotFoundError:
        print(f"O arquivo '{file_path}' não foi encontrado.")
        return False

    return False

# %%
def validate_sequences(file_path: str) -> bool:
    """Verifica se todos os caracteres das sequências são caracteres válidos

    Args:
        file_path (str): Caminho para o arquivo

    Returns:
        bool: True para sequencia válida e False para sequência inválida
    """
    
    valid_characters = set('ACDEFGHIKLMNPQRSTVWY')
    try:
        with open(file_path, 'r') as file:
            for line in file:
                if line.startswith('>'):
                    continue  # Pula a linha de cabeçalho
                sequence = line.strip()
                if not set(sequence).issubset(valid_characters):
                    return False
    except FileNotFoundError:
        print(f"O arquivo '{file_path}' não foi encontrado.")
        return False

    return True

# %%
def remove_pipe(name: str, path_in_fasta: str, input_path: str) -> None:
    """Cria um arquivo de sequências únicas, a partir
    de um arquivo de entrada fasta.

    Args:
        output_file (str): Caminho do arquivo de saída
        path_in_fasta (str): Caminho do arquivo com as sequências de entrada

    Returns:
        str: Caminho do arquivo de saída
    """

    sequences = list(SeqIO.parse(path_in_fasta, "fasta"))
    
    # Criar um dicionário para armazenar as sequências únicas
    unique_sequences = {}

    # Iterar pelas sequências do arquivo de entrada
    for sequence in sequences:
        # Verificar se a sequência já existe no dicionário de sequências únicas
        if str(sequence.seq) not in unique_sequences:
            # Se a sequência é única, armazená-la no dicionário
            unique_sequences[str(sequence.seq)] = sequence

    # Criar uma lista de sequências únicas
    unique_sequences_list = list(unique_sequences.values())

    # Substitui o arquivo original por um arquivo tratado
    SeqIO.write(unique_sequences_list, os.path.join(input_path, name), "fasta")

# %%
def v_sequences(input_path: str):
    """Percorre os arquivos na pasta de entrada.
    Em caso de arquivos com sequências duplicadas ou sequências inválidas, substitui por um arquivo tratado (no_pipe)

    Args:
        path_in_fasta (str): Pasta com os arquivos de entrada em formato fasta
    """

    for name_file in os.listdir(input_path):
        c_path = os.path.join(input_path, name_file)

        if duplicate_names(c_path) or not(validate_sequences(c_path)):
            remove_pipe(Path(name_file), c_path, input_path)

# %%
def clean_files(dir_path: str) -> None:
    """Apaga todas os arquvos de uma pasta
    Obs.: Mantém os arquivos [file.gitkeep e NoPipe]

    [tmp, Trees, full_dataset_plasmodium] - Essas as pastas que normalmente tem que ser limpas

    Args:
        data_output_path (str): caminho da pasta
    """

    files = os.listdir(dir_path)
    if "full_dataset_plasmodium" in dir_path:
        for file in files:
            file_path = os.path.join(dir_path, file)
            if 'NoPipe' in file:
                os.remove(file_path)
    else:       
        for file in files:
            file_path = os.path.join(dir_path, file)
            if file != "file.gitkeep":
                os.remove(file_path)

# %%
def construir_arvores(path_out_aln: str, path_out_tree: str, evolutionary_model:str = 'nj', output_format: str = 'nexus', distance_method: str = 'identity') -> None:
    """_summary_

    Args:
        path_out_aln (str): _description_
        path_out_tree (str): _description_
        evolutionary_model (str, optional): Pode ser "nj" ou "upgm". Defaults to 'nj'.
        output_format (str, optional): Pode ser "newick", "nexus" ou "phyloxml". Defaults to 'nexus'.
        distance_method (str, optional): "identity", "blastn", "trans", "hamming", "kimura", "jukes-cantor", "logdet", "mcc", "poisson", "similarity". Defaults to 'identity'.
    """
    
    clean_files(path_out_tree) # Apaga todos os arquivos de árvores da pasta de saída que estejam lá de execuções anteriores

    for file_aln in os.listdir(path_out_aln):
        print(file_aln)

        if not file_aln.endswith('.aln'):   # Verifica se é um arquivo de alinhamento
            continue
        
        try:
            # Abre o arquivo de alinhamento
            with open(os.path.join(path_out_aln, file_aln), "r") as handle:
                alignment = AlignIO.read(handle, "clustal") # O objeto MultipleSeqAlignment retornado é armazenado na variável.
        except Exception as e:
            print(e)
            continue

        sequence_names = [record.id for record in alignment]
        duplicates = [item for item, count in Counter(sequence_names).items() if count > 1]

        if duplicates:
            print("Nomes duplicados encontrados:", duplicates)

            for i, record in enumerate(alignment):
                record.id = f"seq_{i}"
        
        # Calcula a matriz de distância
        # argumento 'identity', que indica que a distância entre as sequências será medida pelo número de identidades, 
        # ou seja, a fração de posições nas sequências que possuem o mesmo nucleotídeo ou aminoácido.

        calculator = DistanceCalculator(distance_method)

        # Calcula a matriz de distâncias entre as sequências
        distance_matrix = calculator.get_distance(alignment) 

        # Constrói a árvore filogenética
        # Constrói árvores filogenéticas a partir de matrizes de distâncias entre sequências.
        constructor = DistanceTreeConstructor()
        
        match evolutionary_model.lower():
            case 'nj':
                # Para NJ
                tree = constructor.nj(distance_matrix)
            case 'upgma':
                # Para UPGMA
                tree = constructor.upgma(distance_matrix)

        # Salva a árvore
        path_o_tree = os.path.join(path_out_tree,f'tree_{Path(file_aln).stem}.{output_format}')
        Phylo.write(tree, path_o_tree, output_format)

# %%
def align_sequence(
    algoritmo: str,
    path_in_fasta: str, 
    path_out_aln: str, 
    *args, **kwargs
):
    """_summary_

    Args:
        path_in_fasta (str): _description_
        path_out_aln (str): _description_
        path_out_dnd (str): _description_
        path_old_dnd (str): _description_

    Returns:
        _type_: _description_
    """
    input_path, file_name = os.path.split(path_in_fasta)
    file_out_aln = os.path.join(path_out_aln, f'{Path(file_name).stem}.aln')

    match algoritmo.lower():
        case 'muscle':
            command = make_muscle(path_in_fasta, file_out_aln, *args, **kwargs)
            p = subprocess.run(command, capture_output=True, text=True)
        
        case 'clustalw':
            command = make_clustalw(path_in_fasta, file_out_aln, *args, **kwargs)
            p = subprocess.run(command, capture_output=True, text=True)
        
        case 'clustalo':
            # command = make_clustalo(path_in_fasta, file_out_aln, 'auto', 'force', outfmt='clu', threads=4, log='log_file.txt')
            command = make_clustalo(path_in_fasta, file_out_aln, *args, **kwargs)
            p = subprocess.run(command, capture_output=True, text=True)
        
        case 'mafft':
            # command = make_mafft(path_in_fasta, 'auto', retree=2)
            command = make_mafft(path_in_fasta, *args, **kwargs)
            p = subprocess.run(command, capture_output=True)

            with open(file_out_aln, 'wb') as output_file:
                output_file.write(p.stdout)

            # to_clustalw(file_out_aln)
            
        case 'probcons':
            command = make_probcons(path_in_fasta, *args, **kwargs)
            p = subprocess.run(command, capture_output=True)
            
            with open(file_out_aln, 'wb') as output_file:
                output_file.write(p.stdout)

        case 't_coffee':
            command = make_t_coffee(path_in_fasta, file_out_aln, *args, **kwargs)
            p = subprocess.run(command, capture_output=True, text=True)


    if p.stderr:
        print(command)
        print(p.stderr)
        return None

    # Mover o arquivo de saída .dnd para o diretório "resultados"
    file_old_dnd = os.path.join(input_path, f'{Path(file_name).stem}.dnd')
    file_out_dnd = os.path.join(path_out_aln, f'{Path(file_name).stem}.dnd')

    if os.path.exists(file_old_dnd):
        os.rename(file_old_dnd, file_out_dnd)    

# %%
def sub_tree(path: str, name_subtree: str, data_format: str, data_output_path: str,  extension_format: str) -> list:
    """Gera as subarvores a partir de um arquivo de alinhamento .aln

    Args:
        path (str): _description_
        name_subtree (str): _description_
        data_format (str): _description_
        data_output_path (str): _description_
        extension_format (str): _description_

    Returns:
        list: retorna uma lista com todos os arquivos de subarvores do arquivo de entrada
    """
    tree = Phylo.read(path, data_format)
    name_subtree = name_subtree.rsplit(".", 1)[0]

    #Lista caminhos das subárvores (que posteriormente serão utilizadas para compor a matriz de subárvores)
    row_subtree = []

    for clade in tree.find_clades():
        subtree = Phylo.BaseTree.Tree(clade)
        if subtree.count_terminals() > 1:
            filepath_out = os.path.join(data_output_path, f'{name_subtree}_{clade.name}.{extension_format}')
            Phylo.write(subtree, filepath_out, data_format)        
            row_subtree.append(filepath_out)
            
    return row_subtree 

# %%
def directory_has_single_file(directory_path: str) -> str:
    """Verifica se o diretório possui somente um arquivo

    Args:
        directory_path (str): Diretório

    Returns:
        _type_: True para somente um arquivo, False para diferente de um arquivo ou diretório inexistente
    """
    if not os.path.isdir(directory_path):
        return False

    files = os.listdir(directory_path)
    if len(files) != 1:
        return False

    return True

# %%
def preencher_matriz(matriz: list, max_columns:int, valor_preenchimento=None):
    # Preencher as células vazias com o valor de preenchimento
    for row in matriz:
        while len(row) < max_columns:
            row.append(valor_preenchimento)

    return matriz

# %%
def grade_maf(path_1:str, path_2:str, data_format:str) -> int:
    if(path_1 is None or path_2 is None):
        return -1      
    grau = 0

    subtree_1 = Phylo.read(path_1,data_format)
    subtree_2 = Phylo.read(path_2,data_format)

    # Lista todas as clades ( folhas )
    list_1 = [i.name for i in subtree_1.get_terminals()]
    list_2 = [i.name for i in subtree_2.get_terminals()]

    sorted_list1 = sorted(list_1)
    sorted_list2 = sorted(list_2)
    
    for i in range(len(list_1)):
        for j in range(len(list_2)):
            if sorted_list1[i] == sorted_list2[j]:
                grau += 1
    
    return grau

# %%
def fill_dict(dict, max_columns):
    for i in range(max_columns):
        dict[i+1] = {}

    return dict

# %%
def compare_subtrees(max_rows: int, max_columns: int, matrix_subtree: list, dict_maf_database: dict):
    max_maf = 0
    for i in range(max_rows):
        for j in range(max_columns):
            for k in range(max_rows):
                for l in range(max_columns): 
                    if i != k:
                        if max_maf <= grade_maf(matrix_subtree[i][j],matrix_subtree[k][l], 'nexus'):
                            max_maf = grade_maf(matrix_subtree[i][j],matrix_subtree[k][l], 'nexus')

                        g_maf = grade_maf(matrix_subtree[i][j], matrix_subtree[k][l], 'nexus')
                        if g_maf is not False and g_maf >= 1:
                            if g_maf not in dict_maf_database:
                                dict_maf_database[g_maf] = {}
                            if matrix_subtree[i][j] not in dict_maf_database[g_maf]:
                                dict_maf_database[g_maf][matrix_subtree[i][j]] = []
                            dict_maf_database[g_maf][matrix_subtree[i][j]].append(matrix_subtree[k][l])

    return max_maf, dict_maf_database

# %%
def files_align(algoritmo: str, input_path: str, path_out_aln: str, *args, **kwargs):
    clean_files(path_out_aln)

    all_files = os.listdir(input_path)
    files = random.choices(all_files, k=random.random(10, len(all_files)))

    files = os.listdir(input_path)
        
    session = Session()
    for file in files:
        colunas_filtradas = (
            session.query(Entrada.id)
                .filter(Entrada.nome == 'PLASMODIUM6.fasta')
                .first()
            )

        te = Tarefas_Entradas(idTarefa=id_tarefa, idEntrada=colunas_filtradas[0])
        session.add(te)
        session.commit()
        id_tarefa = tarefa.id

        align_sequence(algoritmo, os.path.join(input_path, file), path_out_aln, *args, **kwargs)
    session.close()

# %%
def make_matrix(input_path: str, data_output_path: str, output_format: str) -> tuple:
    clean_files(data_output_path)
    
    files = os.listdir(input_path)

    matrix_subtree = []

    for name_file in files:
        if name_file != "file.gitkeep":
            file_path = os.path.join(input_path, name_file)
            matrix_subtree.append(sub_tree(file_path, name_file, 'nexus', data_output_path, output_format))

    max_columns = max(len(row) for row in matrix_subtree)
    max_rows = len(matrix_subtree)

    # Preencher a matriz
    matrix_subtree = preencher_matriz(matrix_subtree, max_columns, None)

    return matrix_subtree, max_columns, max_rows

# %%
def extrair_informacoes_fasta(input_path: str):
    infos_entradas = []

    files = [file for file in os.listdir(input_path) if os.path.isfile(os.path.join(input_path, file)) and os.path.getsize(os.path.join(input_path, file)) > 1024]
    for arquivo_fasta in files:
        info = {
            "nome_arquivo": arquivo_fasta,
            "numero_sequencias": 0,
            "maior_comprimento": 0,
            "menor_comprimento": float("inf"),
            "comprimento_medio": 0
        }
        
        comprimentos = []
        
        with open(os.path.join(input_path, arquivo_fasta), "r") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                comprimento = len(record.seq)
                comprimentos.append(comprimento)
                info["numero_sequencias"] += 1
                # info["identificadores"].append(record.id)
                info["maior_comprimento"] = max(info["maior_comprimento"], comprimento)
                info["menor_comprimento"] = min(info["menor_comprimento"], comprimento)
        
        # Calcular o comprimento médio
        if info["numero_sequencias"] > 0:
            info["comprimento_medio"] = sum(comprimentos) / info["numero_sequencias"]


        create_or_retrieve(
            Entrada(
                nome=arquivo_fasta, 
                tamanho=os.path.getsize(os.path.join(input_path, arquivo_fasta)),
                qtdSequencias = info['numero_sequencias'],
                maiorComprimento = info['maior_comprimento'],
                menorComprimento = info['menor_comprimento'],
                comprimentoMedio = info['comprimento_medio']
                ),
            Entrada,
            ['nome']
        )

        infos_entradas.append(info)
    
    return infos_entradas

# %%
def salvar_parametros(*args, **kwargs):
    dici = {}
    for arg in args:
        dici[arg] = None
    
    for item, value in kwargs.items():
        dici[item] = value

    return dici

# algoritmos = ['muscle', 'clustalw', 'clustalo', 'mafft', 'probcons', 't_coffee']
algoritmos = ['muscle', 'clustalw']
if __name__ == '__main__':
    for _ in range(300):
        # algoritmo = random.choice(algoritmos)
        algoritmo = 'probcons'
        # %% [markdown]
        # ## 1. Sciphy

        # %% [markdown]
        # ### 1.1 Validação das sequências
        # #### Verifica se no arquivo existem sequências duplicadas e se todos os caracteres da sequência são válidos
        # #### Se houver sequencia duplicadas cria um novo arquivo com sufixo _nopipe e faz as etapas posteriores em cima desse arquivo ao invés do original

        # %%
        v_sequences(os.path.join('data', 'full_dataset_plasmodium'))

        # Coleta informações sobre os arquivos de entrada e já coloca no banco de dados
        infos_entradas = extrair_informacoes_fasta(os.path.join("data", "full_dataset_plasmodium"))

        # %% [markdown]
        # #### Inicia o monitoramento de recursos

        # %%
        # Coleta as informações do Host de execução
        host = create_or_retrieve(
            Host(
                nome=os.uname().nodename, 
                processador=get_cpu_model(),
                capacidade_memoria=psutil.virtual_memory().total / (1024 ** 3)
            ),
            Host,
            ['nome']
        )

        # %%
        initial_disk_io = psutil.disk_io_counters()
        monitor = Execucao(
            HoraInicio=time.time(),
            UsoCPU=psutil.cpu_percent(),
            MemoriaDisponivel=psutil.virtual_memory().free,
        )

        session = Session()
        session.add(monitor)
        session.commit()
        id_monitor = monitor.id
        session.close()

        # %%
        tarefa = Tarefa(nome='Subárvores Frequentes', algoritmo='NMFSt.P', idExecucao=id_monitor, idHost=host.id)
        session = Session()
        session.add(tarefa)
        session.commit()
        id_tarefa = tarefa.id
        session.close()

        # %% [markdown]
        # ### 1.2 Alinhamento múltiplo de sequências

        # %%
        params, tags = sort_params(algoritmo)
                
        d_parametros = salvar_parametros(*tags, **params)
        d_parametros['algoritmo'] = algoritmo

        try:
            # %%
            print(d_parametros)
            files_align(algoritmo, os.path.join('data', 'full_dataset_plasmodium'), os.path.join('data', 'out', 'tmp'), *tags, **params)

            # %%
            session = Session()
            for chave, valor in d_parametros.items():
                parametros = Parametros(Chave=chave, Valor=valor, idTarefa=id_tarefa)
                session.add(parametros)
            session.commit()
            session.close()

            # %% [markdown]
            # ### 1.3 Escolha do modelo evolutivo

            # %% [markdown]
            # 

            # %% [markdown]
            # ### 1.4 Geração da Árvore Filogenética
            # Parametros de entrada do algoritmo: <br>
            # Modelo evoltivo: Pode ser "nj" ou "upgma" <br>
            # Formato de saída: Pode ser nexus <br>

            # %%
            d_parametros = {
                'evolutionary_model':'nj', 
                'output_format':'nexus', 
                'distance_method':'identity'
            }

            session = Session()
            for chave, valor in d_parametros.items():
                parametros = Parametros(Chave=chave, Valor=valor, idTarefa=id_tarefa)
                session.add(parametros)
            session.commit()
            session.close()

            # %%
            print("Construindo árvores: ")
            construir_arvores(os.path.join("data", "out", "tmp"), os.path.join("data", "out", "Trees"), 
                            d_parametros['evolutionary_model'], d_parametros['output_format'], d_parametros['distance_method'])

            # %% [markdown]
            # ### 1.5 Geração das Subárvores Possíveis

            # %%
            matrix_subtree, max_columns, max_rows = make_matrix(os.path.join("data", "out", "Trees"), os.path.join("data", "out", "Subtrees"), "nexus")

            # %% [markdown]
            # ### 1.6 Mapeamento das Subárvores

            # %%
            dict_maf_database = {}
            dict_maf_database = fill_dict(dict_maf_database, max_columns)
                

            # %% [markdown]
            # ### 1.7 Cálculo da Similaridade entre as Subárvores

            # %%
            print("Comparando subárvores")
            max_maf, dict_maf_database = compare_subtrees(max_rows, max_columns, matrix_subtree, dict_maf_database)

            # %% [markdown]
            # ### 1.8 Geração do Dicionário de Saída

            # %%
            print(max_maf)
            for i, j in dict_maf_database.items():
                print(i,j)
                for key, val in j.items():
                    # print(i, key, val)
                    continue

            # %% [markdown]
            # #### Complementando as variáveis de monitoramento

            # %%
            session = Session()
            current_disk_io = psutil.disk_io_counters()
            entrada = session.query(Execucao).filter_by(id=id_monitor).update({
                'LeituraDisco': current_disk_io.read_bytes - initial_disk_io.read_bytes,
                'EscritaDisco': current_disk_io.write_bytes - initial_disk_io.write_bytes,
                'HoraFim': time.time()
            })

            session.commit()
            session.close()

            print("Fim")

        except Exception as e:
            print("ERRO")
            print(e)
            logging.error(e, exc_info=True)