#!/usr/bin/env python
# coding: utf-8

# In[1]:


import sqlite3
import os
import psutil
import shutil
import time
import pandas as pd
import sqlalchemy
import subprocess


# In[2]:


from Bio import AlignIO, Phylo, SeqIO
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor


# In[3]:


from sqlalchemy.orm import declarative_base, sessionmaker
from sqlalchemy import Column, Integer, String, Float, ForeignKey
from sqlalchemy.exc import SQLAlchemyError


# In[4]:


engine = sqlalchemy.create_engine('sqlite:///dados.db')
Base = declarative_base()


# In[5]:


class Entrada(Base):
    __tablename__ = 'Entrada'

    id = Column(Integer, primary_key=True, autoincrement=True)
    nome = Column(String(50), unique=True, nullable=False)
    tamanho = Column(Float)


# In[6]:


class Host(Base):
    __tablename__ = 'Host'
    
    id = Column(Integer, primary_key=True, autoincrement=True)
    nome = Column(String(50), nullable=False, unique=True)
    processador = Column(String(50), nullable=False)
    capacidade_memoria = Column(Float)
    frequencia_memora = Column(Integer)


# In[7]:


class Tarefa(Base):
    __tablename__ = 'Tarefa'

    id = Column(Integer, primary_key=True, autoincrement=True)
    nome = Column(String(50), nullable=False)
    algoritmo = Column(String(50))
    parametros = Column(String(256))


# In[8]:


class Execucao(Base):
    __tablename__ = 'Execucao'

    id = Column(Integer, primary_key=True, autoincrement=True)
    HoraInicio = Column(Float)  
    HoraFim = Column(Float)
    MemoriaDisponivel = Column(Float)
    UsoCPU = Column(Float)
    LeituraDisco = Column(Float)
    EscritaDisco = Column(Float)
    idEntrada = Column(Integer, ForeignKey('Entrada.id')) 
    idHost = Column(Integer, ForeignKey('Host.id'))
    idTarefa = Column(Integer, ForeignKey('Tarefa.id'))


# In[9]:


Base.metadata.create_all(engine) # Cria as tabelas de acordo com as classes


# In[10]:


def cria_arvore_diretorios():
    if not os.path.isdir('files'):
        os.mkdir('files')

    shutil.rmtree(os.path.join('files', 'output'))
    
    if not os.path.isdir(os.path.join('files', 'output')):
        os.mkdir(os.path.join('files', 'output'))
    if not os.path.isdir(os.path.join('files', 'output', 'arvores_filogeneticas')):
        os.mkdir(os.path.join('files', 'output', 'arvores_filogeneticas'))
    if not os.path.isdir(os.path.join('files', 'output', 'sequencias_alinhadas')):
        os.mkdir(os.path.join('files', 'output', 'sequencias_alinhadas'))

    if not os.path.isdir(os.path.join('files', 'input')):
        os.mkdir(os.path.join('files', 'input'))


# In[11]:


def baixar_entradas():
    pass


# In[12]:


def alinhamento(INPUT_PATH,
                OUTPUT_PATH,
                INPUT_FILE, 
                algoritmo):

    OUTPUT_FILE = os.path.join(OUTPUT_PATH, INPUT_FILE.split('.')[0] + '_' + algoritmo + '.aln')
    
    print(f' --- {INPUT_FILE} - {algoritmo} ---')
    match algoritmo.lower():
        case 'clustalw':
            p = subprocess.run(['clustalw', '-INFILE=' + os.path.join(INPUT_PATH, INPUT_FILE), '-OUTFILE=' + OUTPUT_FILE], 
                               capture_output=True)
    
        case 'muscle':
            p = subprocess.run(['muscle', '-align', os.path.join(INPUT_PATH, INPUT_FILE), '-output', OUTPUT_FILE], 
                               capture_output=True)
                
        case 'clustalo':
            p = subprocess.run(['clustalo', '-i', os.path.join(INPUT_PATH, INPUT_FILE), '-o', OUTPUT_FILE],
                               capture_output=True)
    
        case 'mafft':
            with open(OUTPUT_FILE, 'w') as o:
                p = subprocess.run(['mafft', '--auto', os.path.join(INPUT_PATH, INPUT_FILE)], stdout=o)
    
        case 'probcons':
            with open(OUTPUT_FILE, "w") as o:
                p = subprocess.run(['probcons', os.path.join(INPUT_PATH, INPUT_FILE)], stdout=o)
    
        case 't-coffee':
            p = subprocess.run(['t_coffee', '-seq', os.path.join(INPUT_PATH, INPUT_FILE), '-outfile', OUTPUT_FILE], 
                               capture_output=True )

        case other:
            print('Algotitmo não implementado')


# In[13]:


def create_or_retrieve(obj, Classe, atributos):
    session = Session()
        
    filtro = [getattr(Classe, attr) == getattr(obj, attr) for attr in atributos]
    instancia = session.query(Classe).filter(*filtro).first()

    if not instancia:
        session.add(obj)
        session.commit()
        instancia = obj

    while not hasattr(instancia, 'id'):
        time.sleep(1)

    session.expunge(instancia)
    session.close()
    
    return instancia


# In[14]:


def get_cpu_model():
    with open('/proc/cpuinfo') as f:
        for line in f:
            if 'model name' in line:
                return line.split(':')[1].strip()


# In[15]:


algoritmos = ['clustalw', 'muscle', 'clustalo', 'mafft', 'probcons', 't-coffee']
tarefa = 'Sequenciamento'
INPUT_PATH = os.path.join('files', 'input')
OUTPUT_PATH = os.path.join('files', 'output', 'sequencias_alinhadas')
INPUT_FILE = 'ls_orchid.fasta'

Session = sessionmaker(bind=engine)
session = Session()
for INPUT_FILE in os.listdir(INPUT_PATH):
    for algoritmo in algoritmos:
        # Uso de CPU por core. Vai entrar em outra tabela
        # carga_dict = {i: carga for i, carga in enumerate(psutil.cpu_percent(percpu=True))}

        entrada = create_or_retrieve(
            Entrada(nome=INPUT_FILE, tamanho=os.path.getsize(os.path.join(INPUT_PATH, INPUT_FILE))),
            Entrada,
            ['nome']
        )

        host = create_or_retrieve(
            Host(nome=os.uname().nodename, processador=get_cpu_model()),
            Host,
            ['nome']
        )

        tarefa_obj = create_or_retrieve(
            Tarefa(nome=tarefa, algoritmo=algoritmo, parametros=None),
            Tarefa,
            ['nome', 'algoritmo']
        )

        # Pega as variáveis do incio da execução
        initial_disk_io = psutil.disk_io_counters()
        monitor = Execucao(
            HoraInicio=time.time(),
            MemoriaDisponivel=psutil.virtual_memory().free,
            UsoCPU=psutil.cpu_percent(),
            idEntrada=entrada.id,
            idHost = host.id
        )

        # Executa
        alinhamento(INPUT_PATH,
                    OUTPUT_PATH,
                    INPUT_FILE,
                    algoritmo)

        # Pega as variáveis do fim da execução
        current_disk_io = psutil.disk_io_counters()
        monitor.LeituraDisco = current_disk_io.read_bytes - initial_disk_io.read_bytes
        monitor.EscritaDisco = current_disk_io.write_bytes - initial_disk_io.write_bytes
        monitor.HoraFim = time.time()
        monitor.idEntrada = entrada.id
        monitor.idTarefa = tarefa_obj.id
        session.add(monitor)

session.commit()
session.close()


# In[16]:


for file in os.listdir():
    if file.endswith('.dnd'):
        os.remove(file)


# In[ ]:




