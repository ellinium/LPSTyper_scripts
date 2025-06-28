#Genbank to CSV parser
#Set the directory in main()

import os
import csv
from Bio import GenBank
import operator


def gb_2_csv(file_name):

    gb_file_name = file_name
    strain_dir = os.path.split(gb_file_name)[0]
    strain_name_full = os.path.splitext(gb_file_name)[0]
    strain_name = os.path.split(strain_name_full)[1]

    try:
        genbank = open(gb_file_name).read().split('LOCUS  ')  # opens gene bank file and splits by '//\n' to create list of each genes
    except:
        print('Error: File ' + gb_file_name + ' not found')

    contig_num = len(genbank)

    print(file_name)
    print("\nParsing started")
    output_file = strain_dir + '/' + strain_name + '.csv'
    output_file = open(output_file, 'w')  # opening a file to write the ouput
    output = csv.writer(output_file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)

    output.writerow(
        ['Strain', 'Locus', 'Locus_Tag', 'Protein_name', 'CDS_start', 'CDS_stop', 'CDS_complement', 'GeneID', 'AA_seq', 'Protein_length'])

    # going through all the CONTIGS in the list
    for n in range(1, contig_num):  # 0 index is empty

        test = 'LOCUS  ' + genbank[n].lstrip('\n')  # removing new line of from individual genebank files
        query = open('genbank.txt', 'w')  # creating a genbank file to create query gene bank file
        query.write(test)
        query.close()
        parser = GenBank.RecordParser()  # using biopython function for parsing
        record = parser.parse(open('genbank.txt'))
        locus = record.locus
        nt_seq = (record.sequence).strip('\n')  # stores nucleotide sequence
        nm_and_version = (record.version).strip('\n')  # contains nm and nm_version

        # TODOnm = (nm_and_version.split('.')[0]).strip('\n')
        # TODOnm_version = (nm_and_version.split('.')[1]).strip('\n')

        ############################################################################################
        source = record.features[0]  # contains all the fields of source
        try:
            organism = source.qualifiers[0].value.strip('\n') + ':' + source.qualifiers[2].value.strip('\n')
        except:
            organism = ''

        try:
            chrm = (source.qualifiers[3].value).strip('\n')  # stores chromosome number
        except:
            chrm = ''

        try:
            chrm_map = source.qualifiers[4].value.strip('\n')
        except:
            chrm_map = ''

        geneID = ''
        ############################################################################################
        # print(str(record.features))
        if len(record.features) > 1:
            gene = record.features[1]  # contains all the field of gene

        cds = list()
        for c in range(0, len(record.features)):
            if ('CDS' in record.features[c].key or 'tRNA' in record.features[c].key or 'rRNA' in record.features[c].key or 'tmRNA' in record.features[c].key):
                cds.append(record.features[c])
        #	break
        # else:
        #	continue

        # print ("CDS = " + str(len(cds)))
        for cds_num in range(len(cds)):
            cds_cur = cds[cds_num]  # current cds record

            cds_start_stop = (cds_cur.location).strip('\n')  # stores cds start and stop position
            #print('cds start stop = ' + cds_start_stop)
            if 'join' not in cds_start_stop:
                cds_start = (cds_start_stop.split('..')[0]).strip('\n')
                cds_stop = (cds_start_stop.split('..')[1]).strip('\n')
            else:
                cds_start_stop = cds_start_stop.replace('join(', '')
                cds_start_stop = cds_start_stop.replace(')', '')
                cds_start_lst = cds_start_stop.split(',')
                cds_start = ''
                cds_stop = ''
                for c in cds_start_lst:
                    cds_start = cds_start + '||' + (c.split('..')[0]).strip('\n')
                    cds_stop = cds_stop + '||' + (c.split('..')[1]).strip('\n')

            cds_complement = False
            if operator.contains(cds_start_stop, "complement"):
                cds_start = cds_start.replace('complement(', '')
                cds_stop = cds_stop.replace(')', '')
                cds_complement = True

            # creating a empty dictionary to go through the elements in the CDS and update later if present
            cds_dict = {"gene": '', "product": '',
                        "db_xref": '', "translation": '', 'locus_tag': '', "num_aa": '', "gene_synonym": '', "note": ''}

            for n in range(0, len(cds_cur.qualifiers)):  # going through all the elements in the cds
                for key, value in cds_dict.items():  # looping through the dictionary items to see if present in cds
                    #					if ((key in cds_cur.qualifiers[n].key) or (key in cds_cur.qualifiers[n].value)):
                    if (key + '='  in cds_cur.qualifiers[n].key):
                        keys = str(key)  # storing dictionary key
                        cds_dict[keys] = str(cds_cur.qualifiers[n].value)  # updating dictionary key with values
                        break
                    else:
                        continue
            # TODOnp = cds_dict["protein_id"].split('.')[0]+'"'
            # TODOnp_version = '"'+cds_dict["protein_id"].split('.')[1]
            # hgnc=cds_dict["HGNC"]
            # mim=cds_dict["MIM:"]
            geneid = cds_dict["gene"]

            product = cds_dict["product"]
            synonym = cds_dict["gene_synonym"]
            translation = cds_dict["translation"]
            nt_sequence = record.sequence
            locus_tag = cds_dict["locus_tag"]


            num_aa = ''
            if translation != '': num_aa = len(translation)
            note = cds_dict["note"]

            if ("hypothetical protein" in product and len(note) >2):
                product = note


            #Re-annotates hypothetical proteins from 'note'
            gbk_dict = None
            if len(geneid) == 0 and 'hypothetical protein' not in product:
                #parse product value
                product_tmp = product.replace(']', '')
                gene_info = product_tmp.split('[')

                gbk_dict = {gi.split('=')[0]:gi.split('=')[1] for gi in gene_info if len(gi.split('=')) == 2 }
                geneid = gbk_dict.get('gene')
                geneid = str(geneid).replace(" ", "")
                product = gbk_dict.get('protein')


            # if len(hgnc) !=0:
            # 	hgnc = '"'+hgnc.split(':')[2]
            # if len(mim) !=0:
            # 	mim = '"'+mim.split(':')[1]
            # if len(geneid) !=0:
            # geneid = '"'+geneid.split(':')[1]

            # gvalue = name+','+nm+','+nm_version+','+symbol+','+cds_start+','+cds_stop+',' + hgnc +','+\
            # 	mim+','+cds_dict["EC_number"]+','+geneid+ ','+np+','+np_version+','+synonym+','+\
            # 	translation+','+str(num_aa) +','+str(chrm)+','+chrm_map+','+nt_seq+','+organism+'\n'

            ##print("product = " + str(product))
            ##print("cds_start = " + str(cds_start))
            gvalue = [strain_name, locus, locus_tag.replace("\"", ''), str(product).replace("\"", ''), cds_start, cds_stop, str(cds_complement), str(geneid).replace("\"", ''), str(translation).replace("\"", ''), str(num_aa)]
            ##print(gvalue)
            output.writerow(gvalue)
    print("Parsing completed")
    output_file.close()

def get_file_list_by_ext(directory, ext, name_str = ""):
    files = []
    dir_list = os.walk(directory)
    for r, d, f in dir_list:
        for file in f:
            if 'Store' in file:
                a = 1
            file_ext = str.lower(os.path.splitext(file)[1])
            if file_ext != '' and (("*" + file_ext) == ext):
                if (name_str != "" and name_str in file) or name_str == "":
                    files.append(os.path.join(r, file))
    return files

def main():

    gbk_files = get_file_list_by_ext('DirectoryWithGBK', "*.gbk")
    for file in gbk_files:
        gb_2_csv(file)


if __name__ == '__main__':
    main()