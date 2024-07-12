FROM ubuntu:latest

# Atualize os pacotes e instale o Python e o pip
RUN apt-get update && apt-get install -y \
git \
sudo \
nano \
python3 \
python3-pip \
python3-venv \
clustalw \
muscle \
clustalo \
mafft \
probcons \
t-coffee

# Adicione um novo usuário chamado "andre"
RUN useradd -m -s /bin/bash andre && \
echo "andre:012345" | chpasswd && \
usermod -aG sudo andre

# Defina o diretório de trabalho e dá permissão para o usuário
WORKDIR /home/andre/BioPython

# Copia os arquivos
COPY requirements.txt /home/andre/BioPython
RUN chown -R andre:andre /home/andre

# Defina o usuário padrão para o contêiner
USER andre

# Exponha a porta do Jupyter Notebook
EXPOSE 8888

# Cria um ambiente virtual e instala os pacotes
RUN python3 -m venv /home/andre/BioPython/venv \ 
&& /home/andre/BioPython/venv/bin/pip install -r /home/andre/BioPython/requirements.txt

# Cria um volume para a pasta /home/andre
VOLUME /home/andre

# Ativar venv e iniciar o Jupyter Notebook
CMD /bin/bash -c "source /home/andre/BioPython/venv/bin/activate \
                 && jupyter notebook --ip=0.0.0.0 --port=8888 --allow-root"


# docker build -t biopython .
# docker run --name bioContainer -it -p 8888:8888 -v "biopythonV:/home/andre/BioPython" biopython
# docker start BioPyhton
# docker exec -it bioContainer /bin/bash
# docker run --name testeC -it -p 8888:8888 teste