from Bio import SeqIO as fastaReader
from Bio import SeqUtils as statGC
import fnmatch as padrao
import re as regex
import os

#banco de dados (dicionário) formato fasta
taxon_encodeID_db = {}

#lista de taxons a serem analisados
taxons = ["Aotus", "Microcebus", "Otolemur", "Pan", "Callicebus", "Homo", "Nomascus"]

#Funcao: Cria base de dados formato fasta.
#Recebe como parâmetro o diretório dos arquivos .fas
#Retorna dicionário onde cada chave é composta por: TAXON_ENCODEID e sua respectiva sequência nucleotídica
def create_Taxon_EncodeID_db(diretorio):
    
    #caminho absoluto do diretorio onde estao os arquivos .fas
    path = os.path.abspath(diretorio)

    #Laço for para cada arquivo encontrado na lista de conteudo do diretório
    for nome_arquivo_fas in os.listdir(path):

        #Identifica se o arquivo listado está no formato .fas. Se sim, entao o arquivo é tratado no laço for abaixo
        if padrao.fnmatch(nome_arquivo_fas, "*.fas"):

            #cada fasta do arquivo .fas é tratado individualmente
            for fasta in fastaReader.parse(path + "/" + nome_arquivo_fas, "fasta"):    #leitura de arquivo fasta
                
                #O cabeçalho do fasta está presente na lista de taxons?
                if fasta.id in taxons:                          
                    
                    #O encodeID é retornado do próprio nome do arquivo .fas
                    encodeID = get_encodeIDFromENCODE(nome_arquivo_fas)

                    #O cabeçalho(chave do dicionario) do arquivo fasta é formado pelo fasta.id (Aotus, Pan, Homo,...) mais o ENCODEID
                    head = fasta.id + "_" + encodeID
                    
                    #armazena sequencia fasta de acordo com sua chave
                    taxon_encodeID_db[head] = fasta.seq.upper()

    return taxon_encodeID_db

#Função: Recebe como parâmetro o diretório dos arquivos .txt correspondentes aos resultados do RepeatMasker
#As responsabilidades dessa funcação são extrair os dados do RepeatMasker e, com base nesses dados, recuperar as sequencias 
# alvos de LINES respeitando o filtro de corte de tamanho de bases. Essa função faz chamadas à outras funcões, alem de fazer 
# chamada a funcão responsável por escrever em disco os resultados encontrados.
def createDataFromRepeatMasker(diretorio):
    
    #caminho absoluto do diretorio dos arquivos txt do repeatMasker
    path = os.path.abspath(diretorio)                       
    
    #trata cada arquivo encontado no diretório fornecido pelo usuário
    for arquivo_encode in os.listdir(path):
        
        #Se arquivo tratado é do Repeat Masker
        if padrao.fnmatch(arquivo_encode, "Repeat_Encode*"):
            
            #abre arquivo Repeat Encode para leitura
            arquivo = open(path + "/" + arquivo_encode, 'r') 
            
            for linha in arquivo:                            #ler cada linha do arquivo
                isLine = tem_LINE_L1(linha)                  #busca por padrao LINE/L1 em cada linha, retorna TRUE ou False
                
                if isLine:                                  #Se tiver LINE
                    linha = formata_RepeatMasker(linha)     #formata linha do LINE identificado
                    it = linha.split('\t')                  #cria lista onde cada coluna da linha é um elemento definido
                    taxon = it[4]                           #coluna 5 da linha corresponde ao taxon(Aotus, Microcebus etc)
                    
                    #O taxon lido na linha esta presente na lista taxons?
                    if taxon in taxons:
                        encodeID = get_encodeIDFromRepeatMasker(arquivo_encode) 
                        posicao_inicial = int(it[5])                                #posicial inicial da sequência de um LINE
                        posicao_final = int(it[6])                                  #posicao final da sequência de um LINE
                        tamanho_LINE = posicao_final - posicao_inicial              #tamanho de um LINE
                        
                        if tamanho_LINE >= 5000:                                    #filtro de tamanho em pares de base
                            taxon_encodeID = get_taxon_encodeID(taxon, encodeID)    #chave de busca no taxon_encodeID_db
                            
                            if taxon_encodeID_db.__contains__(taxon_encodeID):         #esse taxon_encodeID está no banco taxon_encodeID_db?
                                
                                #Recupera sequência nucleotídica de acordo com o ID
                                sequencia = get_Seq_From_TaxonEncodeIDdb(taxon_encodeID)

                                #Recupera a sequência nuclotídica do LINE
                                sequencia_LINE = get_seq_LINE(sequencia, posicao_inicial, posicao_final)             
                                
                                #Formata cabecalho do LINE
                                cabecalho_LINE = format_cabecalho_LINE(taxon, encodeID, posicao_inicial, posicao_final) 
                                
                                #formata um arquivo fasta
                                alvo_fasta_LINE = format_fasta_LINE(cabecalho_LINE, sequencia_LINE) 
                                
                                #salva em disco 
                                escreverArquivoFasta(cabecalho_LINE, alvo_fasta_LINE)
                        
            arquivo.close()

#Função: Recebe o caminho relativo do arquivo .fas e retorna apenas nome do arquivo sem a extensao:
# Ex.: fas/ENCODE114.fas --> ENCODE114 
def get_encodeIDFromENCODE(nome_arquivo_fas):
    encodeID = os.path.basename(nome_arquivo_fas).strip(".fas")
    return encodeID

#Função: Recebe um taxon(Aotus, ou Microcebus, ou etc) e encodeID(ENCODE114, ou ENCODE110, ou etc).
#Retorna String: ex.: Aotus_ENCODE110
def get_taxon_encodeID(taxon, encodeID):
    taxon_encodeID = "{}_{}".format(taxon, encodeID)
    return taxon_encodeID

#Função: Recebe uma String, taxon_encodeID, ex.: Aotus_ENCODE110
#Retorna String: sequencia fasta correspondente ao taxon_encodeID
def get_Seq_From_TaxonEncodeIDdb(taxon_encodeID):
    sequencia = taxon_encodeID_db[taxon_encodeID]
    return sequencia

#Funcao: recebe duas Strings, cabecalho de uma arquivo fasta e sua sequencia nucleotídica
#Retona String: formato fasta 
def format_fasta_LINE(cabecalho_LINE, sequencia_LINE):
    fasta_LINE = ">{}\n{}".format(cabecalho_LINE, sequencia_LINE)
    return fasta_LINE

#Funcao. Recebe 4 parâmetros(String, String, int, int) relacionados ao taxon(Aotus, Pan, Homo ou etc), ao ENCODEID(ENCODE110 ou etc),
#e as posicões de inicio e fim em uma sequencia nucleotídica
#Retorna String: cabecalho sequencia fasta formatado. Ex.: (Aotus, ENCODE114, 100, 300) --> Aotus_ENCODE114_100_300
def format_cabecalho_LINE(taxon, encode_id, posicao_inicial, posicao_final):
    cabecalho_fasta_LINE = "{}_{}_{}_{}".format(taxon, encode_id, posicao_inicial, posicao_final)
    return cabecalho_fasta_LINE

#Funcao: Recebe 3 parâmetros(String, int, int) que são uma sequência nucleotídica, 
#posicao inicial e final referentes a localização da sequência nucleotídica do LINE
#Retorna uma substring(sequência nucloetídica) da posicao inicial e final da 
#sequencia nueclotidica de um taxon no banco de dados de fasta
def get_seq_LINE(sequencia_fasta, posicao_inicial, posicao_final):
    sequencia = sequencia_fasta[posicao_inicial:posicao_final]
    return sequencia

#Funcao: Recebe nome do arquivo RepeatMasker e retona apenas o encode_id. Ex.: Repeat_Encode_010.txt para ENCODE010 
def get_encodeIDFromRepeatMasker(nome_arquivo):
    encode_id = regex.sub("Repeat_Encode_", "ENCODE", nome_arquivo).strip(".txt")
    return encode_id

 #Funcao: Retorna TRUE se a linha tiver o padrao LINE/L1, caso contrario, retorna FALSE
def tem_LINE_L1(linha_do_repeat_masker):
    return regex.search("LINE/L1", linha_do_repeat_masker)

#Funcao: Substitui todos os espaços entre as colunas e por unico TAB. Retorna uma linha com elementos separadas por TAB
def formata_RepeatMasker(linha_do_repeat_masker):
    linha_formatada = linha_do_repeat_masker.strip("^[ ]+") #retira todos os espaços antes da 1a coluna 
    linha_formatada = regex.sub('[ ]+', '\t', linha_formatada) #retira todos os espaços entre as colunas e substituiu por TAB
    return linha_formatada

#Funcao: Recebe 2 parâmetros(String, String)
def escreverArquivoFasta(cabecalho_fasta, sequencia_nucleotidica):
    
    filename = cabecalho_fasta + ".fasta"                       #O nome do arquivo será o mesmo do cabeçalho acresentando a extensao .fasta
    arquivo_saida = open(filename, 'w')
    arquivo_saida.write(str(sequencia_nucleotidica))
    arquivo_saida.close()
    print("O arquivo "+ filename + " foi criado!")

#criacao dos diretorios para cada taxon e mover os resultados para os respectivos diretórios    
def make():
    
    print("\n... Organizando os resultados...")
    os.system("mkdir -p Aotus Callicebus Homo Pan Otolemur Nomascus Microcebus")
    print("\nDiretórios criados!")
    print("...Movendo resultados para seus respectivos diretórios.")
    os.system("mv Aotus_ENCODE* Aotus; mv Callicebus_ENCODE* Callicebus; mv Homo_ENCODE* Homo; mv Pan_ENCODE* Pan; mv Otolemur_ENCODE* Otolemur; mv Nomascus_ENCODE* Nomascus; mv Microcebus_ENCODE* Microcebus")   
    print("\nTrabalho finalizado!")

if __name__ == '__main__':
    #execucao: python3 encodes.py <diretorio_dos_arquivos .fas> <diretorio dos arquivos do RepeatMasker>
    #exe.: python3 encodes.py fas/ Repeat_Encode/
    import sys
    import time

    diretorio_encode = sys.argv[1]  #diretório dos arquivos .fas do ENCODE
    diretorio_repeat = sys.argv[2]  #diretório dos arquivos .txt do RepeatMasker
    
    print("..... Criando taxon_encodeID_db .....")
    create_Taxon_EncodeID_db(diretorio_encode)
    print("taxon_encodeID_db criado!")
    print(".... Extraindo LINES ....")
    time.sleep(2)
    createDataFromRepeatMasker(diretorio_repeat)
    print("\nExtração concluída!")
    make()