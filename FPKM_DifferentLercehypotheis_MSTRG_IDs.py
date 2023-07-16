import re
with open("Lerer_commonprots", 'r') as Lerer_noorth_present, \
    open("Lerernofilt_common_lerceproteins_inallfiles", 'r') as Lerernofilt, open("Lerce_gene_gcenexp_matrix.txt", 'r') as Lerce_genecount, \
        open("GFF_File_edited", 'r') as GFF, open("common_lerceproteins_inallfiles", 'r') as Lerce_nolerer, \
        open("Lerce_Nolerer_FPKM", 'w') as NoLererOutput, open("Lerce_Lerer_lessthanagaric_FPKM", 'w') as SomeLererOutput, \
        open("Lerce_lerernofilt_FPKM", 'w') as LerernofiltOutput, \
        open("input_genome.merged.gtf", 'r') as Input_genome:
    next(Lerce_nolerer)
    next(Lerer_noorth_present)
    next(Lerernofilt)

    Lerce_nolerer_set=set()
    gff_file_numbers_nolerer = set()
    Lerce_somelerer_set=set()
    gff_file_numbers_somelerer = set()
    Lerce_Lerernofilt_set=set()
    gff_file_numbers_lerernofilt = set()



    for Ler in Lerce_nolerer:
        Lerd=Ler.strip().split("|", 1)[1]
        Lerce_nolerer_set.add(Lerd)
    for Lersome in Lerer_noorth_present:
        rd = Lersome.strip().split("|", 1)[1]
        Lerce_somelerer_set.add(rd)
    for lerall in Lerernofilt:
        Lce=lerall.strip().split("|", 1)[1]
        Lerce_Lerernofilt_set.add(Lce)
    data_dict = {}
    for Hopefully in Input_genome:
        if not Hopefully.startswith("tig"):
            continue
        mRNAfield = Hopefully.strip().split('gene_id', 1)[1]
        mRNAfield = mRNAfield.strip().split(';')[1]
        if "mRNA" in mRNAfield:
            mRNAfield_putin=mRNAfield.strip().split("transcript_id")[1]
            MSTRG_field = Hopefully.strip().split('gene_id', 1)[1]
            MSTRG_field=MSTRG_field.strip().split(";", 1)[0]
            proteinname=Hopefully.strip().split("|", 1)[1]
            proteinname=proteinname.strip().split(";")[0]
            Value_tuple = (mRNAfield_putin, proteinname)

            if MSTRG_field in data_dict:
                # Check if value tuple already exists in the list of values
                if Value_tuple not in data_dict[MSTRG_field]:
                    data_dict[MSTRG_field].append(Value_tuple)
            else:
                # If key does not exist, create a new list with the value tuple
                data_dict[MSTRG_field] = [Value_tuple]
    quotes_gone_dict={}
    for key, values in data_dict.items():
        cleankey=key.strip('"')
        cleanvalues = [(val[0].strip().strip('"'), val[1].strip().strip('"')) for val in values]
        quotes_gone_dict[cleankey]=cleanvalues

    Lerce_lernotpresent_dict={}
    for Ke, Va in quotes_gone_dict.items():
        Nolerer_filtered_values = [(mRNA, protein) for mRNA, protein in Va if protein in Lerce_nolerer_set]
        if Nolerer_filtered_values:
            Lerce_lernotpresent_dict[Ke] = Nolerer_filtered_values
    print(Lerce_lernotpresent_dict)

    Lerce_lerersomepres_dict={}
    for k, v in quotes_gone_dict.items():
        someLerer_filt=[(mRNA, protein) for mRNA, protein in v if protein in Lerce_somelerer_set]
        if someLerer_filt:
            Lerce_lerersomepres_dict[k]=someLerer_filt
    Lerce_nofiltlerer_dict = {}
    for ky, vs in quotes_gone_dict.items():
        NofiltLerer_filt = [(mRNA, protein) for mRNA, protein in vs if protein in Lerce_Lerernofilt_set]
        if NofiltLerer_filt:
            Lerce_nofiltlerer_dict[ky] = NofiltLerer_filt
    ###now everything is set up, time to filter, there are some doubleups so need to account for that
    header = Lerce_genecount.readline().strip()
    Lerce_genecount_lines = Lerce_genecount.readlines()
    NoLererOutput.write("ProteinID\tmRNAID\t" + header + "\n")
    SomeLererOutput.write("ProteinID\tmRNAID\t" + header + "\n")
    LerernofiltOutput.write("ProteinID\tmRNAID\t" + header + "\n")
    for noler in Lerce_genecount_lines[1:]:
        nolercols = noler.strip().split('\t')
        NoLererGeneID = nolercols[0]
        if NoLererGeneID in Lerce_lernotpresent_dict:
            nolerervalues = Lerce_lernotpresent_dict[NoLererGeneID]
            NoLererprotein_ids = ";".join([value[1] for value in nolerervalues])
            NoLerermrna_ids = ";".join([value[0] for value in nolerervalues])

            NoLererline = "\t".join([NoLererprotein_ids, NoLerermrna_ids] + nolercols) + "\n"
            NoLererOutput.write(NoLererline)
    for someLer in Lerce_genecount_lines[1:]:
        Somelercols = someLer.strip().split('\t')
        SomeLererGeneID = Somelercols[0]
        if SomeLererGeneID in Lerce_lerersomepres_dict:
            somelerervalues = Lerce_lerersomepres_dict[SomeLererGeneID]
            SomeLererprotein_ids = ";".join([value[1] for value in somelerervalues])
            SomeLerermrna_ids = ";".join([value[0] for value in somelerervalues])
            SomeLererline = "\t".join([SomeLererprotein_ids, SomeLerermrna_ids] + Somelercols) + "\n"
            SomeLererOutput.write(SomeLererline)
    for nofiltLer in Lerce_genecount_lines[1:]:
        nofiltlercols = nofiltLer.strip().split('\t')
        nofiltlerGeneID = nofiltlercols[0]
        if nofiltlerGeneID in Lerce_nofiltlerer_dict:
            nofiltlervalues = Lerce_nofiltlerer_dict[nofiltlerGeneID]
            nofiltlerprotein_ids = ";".join([value[1] for value in nofiltlervalues])
            nofiltlermrna_ids = ";".join([value[0] for value in nofiltlervalues])
            nofiltlerline = "\t".join([nofiltlerprotein_ids, nofiltlermrna_ids] + nofiltlercols) + "\n"
            LerernofiltOutput.write(nofiltlerline)

# Open the files for reading and writing
with open("Lerce_Nolerer_FPKM", 'r+') as Lerce_Nolerer_FPKM_removinghashtag, \
     open("Lerce_lerernofilt_FPKM", "r+") as Lerce_lerernofilt_FPKM_hashtag, \
     open("Lerce_Lerer_lessthanagaric_FPKM", 'r+') as Lerce_Lerer_lessthanagaric_FPKM_hashhtag:

    # Read the first lines
    Lerce_Nolerer_FPKM_removinghashtag_firstline = Lerce_Nolerer_FPKM_removinghashtag.readline()
    Lerce_lerernofilt_FPKM_hashtag_firstline = Lerce_lerernofilt_FPKM_hashtag.readline()
    Lerce_Lerer_lessthanagaric_FPKM_hashhtag_firstline = Lerce_Lerer_lessthanagaric_FPKM_hashhtag.readline()

    # Modify the first lines
    MLerce_Nolerer_FPKM_removinghashtag_firstline = Lerce_Nolerer_FPKM_removinghashtag_firstline.replace("#", '')
    MLerce_lerernofilt_FPKM_hashtag_firstline = Lerce_lerernofilt_FPKM_hashtag_firstline.replace("#", '')
    MLerce_Lerer_lessthanagaric_FPKM_hashhtag_firstline = Lerce_Lerer_lessthanagaric_FPKM_hashhtag_firstline.replace("#", "")

    # Move file pointers back to the beginning
    Lerce_Nolerer_FPKM_removinghashtag.seek(0)
    Lerce_lerernofilt_FPKM_hashtag.seek(0)
    Lerce_Lerer_lessthanagaric_FPKM_hashhtag.seek(0)

    # Write the modified lines back to the files
    Lerce_Nolerer_FPKM_removinghashtag.write(MLerce_Nolerer_FPKM_removinghashtag_firstline)
    Lerce_lerernofilt_FPKM_hashtag.write(MLerce_lerernofilt_FPKM_hashtag_firstline)
    Lerce_Lerer_lessthanagaric_FPKM_hashhtag.write(MLerce_Lerer_lessthanagaric_FPKM_hashhtag_firstline)


#now to remove ay