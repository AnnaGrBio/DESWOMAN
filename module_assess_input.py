def openFile(NameFile):
    F=open(NameFile, "r")
    L=F.readlines()
    return L 


def fill_default_values(dico_elts):
    new_dico_elts = {}
    #print ("lili")
    #print (dico_elts["TPM_threeshold"])
    if "TPM_threeshold" in dico_elts.keys():
        if dico_elts["TPM_threeshold"] == 'Default':
            new_dico_elts["TPM_threeshold"] = 0.5
        else:
            new_dico_elts["TPM_threeshold"] = int(dico_elts["TPM_threeshold"])

    if "transcript_overlap" in dico_elts.keys():
        if dico_elts["transcript_overlap"] == "Default":
            new_dico_elts["transcript_overlap"] = ["intergenic"]
        else:
            list_overlap = dico_elts["transcript_overlap"].split(",")
            new_dico_elts["transcript_overlap"] = list_overlap
            #print (list_overlap)
    if "ORFs_definition" in dico_elts.keys():
        if dico_elts["ORFs_definition"] == "Default":
            new_dico_elts["ORFs_definition"] = "forward"

    if "ORFs_choice" in dico_elts.keys():
        if dico_elts["ORFs_choice"] == "Default":
            new_dico_elts["ORFs_choice"] = [["longest"],["duplicate_handle"]]
        else:
            final_list = [["duplicate_handle"]]
            list_elts_orf_choice = dico_elts["ORFs_choice"].split(";")
            for elt in list_elts_orf_choice:
                sublist = elt.split(",")
                final_list.append(sublist)
            new_dico_elts["ORFs_choice"] = final_list
        #print (new_dico_elts["ORFs_choice"])
    if "parameters_database_prot" in dico_elts.keys():
        if dico_elts["parameters_database_prot"] == "Default":
            new_dico_elts["parameters_database_prot"] = {"type" : "blastp", "mode" : "--more-sensitive"}
        else:
            dico_params_prot = {"type" : "blastp"}
            mode = dico_elts["parameters_database_prot"]
            if mode == "--sensitive" or mode == "-more-sensitive" or mode == "--very-sensitive" or mode == "--ultra-sensitive":
                dico_params_prot["mode"] = mode
            else:
                dico_params_prot["mode"] = "--more-sensitive"
            new_dico_elts["parameters_database_prot"] = dico_params_prot
    else:
        new_dico_elts["parameters_database_prot"] = {"type" : "blastp", "mode" : "--more-sensitive"}
    if "parameters_database_nucl" in dico_elts.keys():
        if dico_elts["parameters_database_nucl"] == "Default":
            new_dico_elts["parameters_database_nucl"] = {"type" : "blastn", "e_value" : "0.01", "coverage" : None, "strand" : None}
        else:
            dico_params_nucl = {"type" : "blastn"}
            list_elts = dico_elts["parameters_database_nucl"].split(",")
            if list_elts[0] == "Default":
                dico_params_nucl["e_value"] = "0.01"
            else:
                dico_params_nucl["e_value"] = list_elts[0]
            if list_elts[1] == "Default":
                dico_params_nucl["coverage"] = None
            else:
                dico_params_nucl["coverage"] = int(list_elts[1])
            if list_elts[2] == "Default":
                dico_params_nucl["strand"] = None
            else:
                dico_params_nucl["strand"] = list_elts[2]
            new_dico_elts["parameters_database_nucl"] = dico_params_nucl
    else:
        new_dico_elts["parameters_database_nucl"] = {"type" : "blastn", "e_value" : "0.01", "coverage" : None, "strand" : None}
    if "filter_genic" in dico_elts.keys():
        if dico_elts["filter_genic"] == "Default":
            dico_elts["filter_genic"] = False
        else:
            dico_elts["filter_genic"] = True


    if "synteny_window" in dico_elts.keys():
        if dico_elts["synteny_window"] == "Default":
            new_dico_elts["synteny_window"] = 2
        else:
            new_dico_elts["synteny_window"] = int(dico_elts["synteny_window"])

    if "premature_stop" in dico_elts.keys():
        if dico_elts["premature_stop"] == "Default":
            new_dico_elts["premature_stop"] = 50
        else:
            new_dico_elts["premature_stop"] = int(dico_elts["premature_stop"])

    if "filter_TE" in dico_elts.keys():
        if dico_elts["filter_TE"] == "Default":
            new_dico_elts["filter_TE"] = "False"

    for variable in dico_elts.keys():
        if variable not in new_dico_elts.keys():
            new_dico_elts[variable] = dico_elts[variable]
    return new_dico_elts


def withdraw_strategy_data(strategy_file):
    dico_elts = {}
    for line in strategy_file:
        if line[0] != "#" and "=" in line:
            line = line.split("\n")[0]
            elts_line = line.split(" = ")
            if len(elts_line) == 2:
                dico_elts[elts_line[0]] = elts_line[1]
    new_dico_elts = fill_default_values(dico_elts)
    return new_dico_elts



