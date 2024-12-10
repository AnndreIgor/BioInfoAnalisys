from sqlalchemy.orm import declarative_base, sessionmaker
from sqlalchemy import Column, Integer, String, Float, ForeignKey, create_engine
from sqlalchemy.exc import SQLAlchemyError
import time

engine = create_engine('sqlite:///dados.db')
Base = declarative_base()
Session = sessionmaker(bind=engine)

class Entrada(Base):
    __tablename__ = 'Entrada'

    id = Column(Integer, primary_key=True, autoincrement=True)
    nome = Column(String(50), unique=True, nullable=False)
    tamanho = Column(Float)
    qtdSequencias = Column(Integer)
    maiorComprimento = Column(Integer)
    menorComprimento = Column(Integer)
    comprimentoMedio = Column(Float)

class Host(Base):
    __tablename__ = 'Host'
    
    id = Column(Integer, primary_key=True, autoincrement=True)
    nome = Column(String(50), nullable=False, unique=True)
    processador = Column(String(50), nullable=False)
    capacidade_memoria = Column(Float)
    frequencia_memora = Column(Integer)
    teraflops = Column(Float)

class Tarefa(Base):
    __tablename__ = 'Tarefa'

    id = Column(Integer, primary_key=True, autoincrement=True)
    nome = Column(String(50), nullable=False)
    algoritmo = Column(String(50))
    idExecucao = Column(Integer, ForeignKey('Execucao.id'))
    idHost = Column(Integer, ForeignKey('Host.id')) 

class Execucao(Base):
    __tablename__ = 'Execucao'

    id = Column(Integer, primary_key=True, autoincrement=True)
    MemoriaDisponivel = Column(Float)
    UsoCPU = Column(Float)
    LeituraDisco = Column(Float)
    EscritaDisco = Column(Float)
    HoraInicio = Column(Float)  
    HoraFim = Column(Float)

class Parametros(Base):
    __tablename__ = 'Parametros'

    id = Column(Integer, primary_key=True, autoincrement=True)
    Chave = Column(String(30), nullable=False)  
    Valor = Column(String(50))
    idTarefa = Column(Integer, ForeignKey('Tarefa.id'))

class Tarefas_Entradas():
    __tablename__ = 'Tarefas_Entradas'

    idTarefa = Column(Integer, ForeignKey('Tarefa.id'))
    idEntrada = Column(Integer, ForeignKey('Entrada.id'))

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

if __name__ == '__main__':
    print("Criando Tabelas")
    Base.metadata.create_all(engine) # Cria as tabelas de acordo com as classes
