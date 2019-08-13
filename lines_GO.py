from Bio import SeqIO as fastaReader
from Bio import SeqUtils as calculate
import fnmatch as padrao
import re as regex
import os
from operator import itemgetter

#banco de dados (dicionário) formato fasta
taxon_encodeID_db = {}
taxonStat = {}
freq_taxon_encode = {}

#lista de taxons a serem analisados
taxons = ["Aotus", "Microcebus", "Otolemur", "Pan", "Callicebus", "Homo", "Nomascus", "Callithrix"]

#Funcao: Cria base de dados formato fasta.
#Recebe como parâmetro o diretório dos arquivos .fas
#Retorna dicionário onde cada chave é composta por: TAXON_ENCODEID e sua respectiva sequência nucleotídica
def create_Taxon_EncodeID_db(diretorio):
    
    path = os.path.abspath(diretorio)

    for nome_arquivo_fas in os.listdir(path): 
        
        if padrao.fnmatch(nome_arquivo_fas, "*.fas"):   #Identifica se o arquivo listado está no formato

            for fasta in fastaReader.parse(path + "/" + nome_arquivo_fas, "fasta"): #para cada fasta do arquivo .fas
                
                if fasta.id in taxons: 
                    
                    encodeID = get_encodeIDFromENCODE(nome_arquivo_fas)
                    head = fasta.id + "_" + encodeID
                    taxon_encodeID_db[head] = fasta.seq.upper() #armazena sequencia fasta de acordo com sua chave

    return taxon_encodeID_db

#Função: Recebe como parâmetro o diretório dos arquivos .txt correspondentes aos resultados do RepeatMasker
#As responsabilidades dessa funcação são extrair os dados do RepeatMasker e, com base nesses dados, recuperar as sequencias 
# alvos de LINES respeitando o filtro de corte de tamanho de bases. Essa função faz chamadas à outras funcões, alem de fazer 
# chamada a funcão responsável por escrever em disco os resultados encontrados.
def createDataFromRepeatMasker(diretorio):
    
    
    path = os.path.abspath(diretorio)                       
    
    
    for arquivo_encode in os.listdir(path):
        
        
        if padrao.fnmatch(arquivo_encode, "Repeat_Encode*"):  #Se arquivo tratado é do Repeat Masker
            
            arquivo = open(path + "/" + arquivo_encode, 'r') 
            
            for linha in arquivo:                            #ler cada linha do arquivo
                isLine = tem_LINE_L1(linha)                  #busca por padrao LINE/L1 em cada linha, retorna TRUE/FALSE
                
                if isLine:                                  #Se tiver LINE
                    linha = formata_RepeatMasker(linha)     #formata linha do LINE identificado
                    it = linha.split('\t')                  #cria lista onde cada coluna da linha é um elemento definido
                    taxon = it[4]                           #coluna 5 da linha corresponde ao taxon(Aotus,Pan e etc)
                    
                    if taxon in taxons:                                             #O taxon lido está presente na lista
                        encodeID = get_encodeIDFromRepeatMasker(arquivo_encode) 
                        posicao_inicial = int(it[5])                                #posicial inicial seq de um LINE
                        posicao_final = int(it[6])                                  #posicao final da seq de um LINE
                        tamanho_LINE = posicao_final - posicao_inicial              #tamanho de um LINE
                        
                        if tamanho_LINE >= 5000:                                    #filtro de tamanho em pares de base
                            taxon_encodeID = get_taxon_encodeID(taxon, encodeID)    #chave de busca no taxon_encodeID_db
                            
                            if taxon_encodeID_db.__contains__(taxon_encodeID):      #taxon_encodeID está no banco? 
                                
                                sequencia = get_Seq_From_TaxonEncodeIDdb(taxon_encodeID)
                                sequencia_LINE = get_seq_LINE(sequencia, posicao_inicial, posicao_final)
                                cabecalho_LINE = format_cabecalho_LINE(taxon, encodeID, posicao_inicial, posicao_final)
                                alvo_fasta_LINE = format_fasta_LINE(cabecalho_LINE, sequencia_LINE)
                                
                                escreverArquivoFasta(cabecalho_LINE, alvo_fasta_LINE)
                                taxonStat[cabecalho_LINE] = calculate.GC(alvo_fasta_LINE)  
                                updateFreqTaxonEncode(taxon, encodeID)                      
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
    
    filename = cabecalho_fasta + ".fasta"              #O nome do arquivo será o mesmo do cabeçalho acresentando a extensao .fasta
    arquivo_saida = open(filename, 'w')
    arquivo_saida.write(str(sequencia_nucleotidica))
    arquivo_saida.close()
    print("O arquivo "+ filename + " foi criado!")


def updateFreqTaxonEncode(taxon, encode):
    cabecalho = taxon + "," + encode

    if freq_taxon_encode.__contains__(cabecalho):
        freq_taxon_encode[cabecalho] += 1
    else:
        freq_taxon_encode[cabecalho] = 1

def relatorioGC():
    report = "SEQ,TAMANHO_SEQ,GC\n"
        
    for head, gc in taxonStat.items():
        aux = head.split("_")
        tamanho = int(aux[3])-int(aux[2]) + 1
        report += "{},{},{:.2f}\n".format(head, tamanho, gc)
    return report           

def relatorioFreq():
    report = "TAXON,ENCODE,FREQUENCIA\n"
    sorted_freq = sorted(freq_taxon_encode.items(), key = itemgetter(1), reverse = True)

    for k,v in sorted_freq:
        report += "{},{}\n".format(k,v)
    return report

def salvarRelatorioFreq(report):
    filename = "taxon_encode_count.csv"
    writer = open(filename, 'w')
    writer.write(str(report))
    writer.close()
    print("O arquivo " + filename + "foi gerado")

def salvarRelatorioGC(report):
    filename = "reportGC.csv"
    writer = open(filename, 'w')
    writer.write(str(report))
    writer.close()
    print("O arquivo " + filename + "foi gerado")

#criacao dos diretorios para cada taxon e mover os resultados para os respectivos diretórios    
def make():
    
    print("\n... Organizando os resultados...")
    os.system("mkdir -p Aotus Callicebus Homo Pan Otolemur Nomascus Microcebus Callithrix")
    print("\nDiretórios criados!")
    print("...Movendo resultados para seus respectivos diretórios.")
    os.system("mv Aotus_ENCODE* Aotus; mv Callicebus_ENCODE* Callicebus; mv Homo_ENCODE* Homo; mv Pan_ENCODE* Pan; mv Otolemur_ENCODE* Otolemur; mv Nomascus_ENCODE* Nomascus; mv Microcebus_ENCODE* Microcebus; mv Callithrix_ENCODE* Callithrix")   
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
    salvarRelatorioGC(relatorioGC())
    salvarRelatorioFreq(relatorioFreq())