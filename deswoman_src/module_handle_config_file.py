import re
import os
from module_colors import *


__author__ = "Anna Grandchamp"
__contributor__=""
__copyright__ = ""
__credits__ = []
__license__ = ""
__version__="1.0.0"
__maintainer__ = "Anna Grandchamp"
__email__ = "anna.grandchamp@inserm.fr"


def parse_config(file_path):
    """
    Parses a configuration file and extracts key-value pairs.

    This function reads a configuration file and extracts key-value pairs from it. It supports different 
    types of values including integers, lists, strings, and `None` for the keyword "Default". The function 
    also handles inline comments and ignores empty lines or lines starting with `#`.

    Args:
        file_path (str): The path to the configuration file to be parsed.

    Returns:
        dict: A dictionary where the keys are the configuration parameters (as strings) and the values 
              are the parsed values (integers, strings, lists, or `None`).

    Example:
        If the configuration file `config.txt` contains:
            ```
            # This is a comment
            param1 = 10
            param2 = "string_value"
            param3 = "1, 2, 3"
            param4 = Default
            ```

        The function will return:
            {
                "param1": 10,
                "param2": "string_value",
                "param3": ["1", "2", "3"],
                "param4": None
            }
    """
    config = {}
    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            # Skip empty lines and comments
            if not line or line.startswith("#"):
                continue  
            # Remove inline comments (// ...)
            line = re.sub(r'//.*', '', line).strip()
            # Extract key-value pairs
            match = re.match(r'(\w+)\s*=\s*(.+)', line)
            if match:
                key, value = match.groups()
                value = value.strip()
                # Convert integers
                if value.isdigit():
                    value = int(value)
                # Convert lists (comma-separated values)
                elif "," in value:
                    value = [v.strip() for v in value.split(",")]
                # Convert "Default" into a Python-friendly value
                elif value.lower() == "default":
                    value = None
                # Convert unquoted strings (file paths, gene names)
                elif not (value.startswith('"') and value.endswith('"')):
                    value = value.strip()
                else:
                    value = value.strip('"')  # Remove surrounding quotes
                config[key] = value
    return config


def string_checking(value_param : str) -> bool:
    """
    Checks if the provided parameter is a string.

    This function verifies whether the given `value_param` is of type string. If it is, the function 
    returns `True`, indicating that the parameter is valid. If the parameter is not a string, it 
    prints an error message and returns `False`.

    Args:
        value_param (str): The value to check, expected to be a string.

    Returns:
        bool: `True` if the value is a string, `False` otherwise.
    """
    if type(value_param) == str:
        return True
    else:
        print (BRIGHT_RED + "ERROR... THE PARAMETER " + RESET + str(value_param) + BRIGHT_RED + " MUST BE A STRING ..." + RESET)
        return False


def check_true_repository(dico_variables : dict, param : str, dico_config : dict) -> bool:
    """
    Validates whether a directory exists at the specified path in the configuration.

    This function checks whether the path provided for a given parameter in the `dico_config` dictionary 
    points to an existing directory. If the directory exists, it updates the `dico_variables` dictionary 
    with the directory path (removes any trailing slashes). If the directory does not exist or the path 
    is invalid, it returns `False` and prints an error message.

    Args:
        dico_variables (dict): A dictionary to store valid directory paths for parameters.
        param (str): The key in `dico_config` corresponding to the directory path to validate.
        dico_config (dict): A dictionary containing configuration values, where each key corresponds 
                             to a parameter and its associated value.

    Returns:
        bool: `True` if the directory exists and is valid, `False` otherwise.
    """
    link = dico_config[param]
    valid_parameter = string_checking(link)
    if valid_parameter == True:
        if os.path.isdir(link) == False:
            print (BRIGHT_RED + "ERROR... THE FOLDER " + RESET + link + BRIGHT_RED + " DOES NOT EXIST ..." + RESET)
            return False
        else:
            if link[len(link) -1] == "/":
                link = link[0:len(link) -1]
            dico_variables[param] = link
            return True
    else:
        return False


def check_true_file(dico_variables : dict, param : str, dico_config : dict) -> bool:
    """
    Validates whether a file exists at the specified path in the configuration.

    This function checks whether the path provided for a given parameter in the `dico_config` dictionary 
    points to an existing file. If the file exists, it updates the `dico_variables` dictionary with the 
    file path. If the file does not exist or the path is invalid, it returns `False` and prints an error message.

    Args:
        dico_variables (dict): A dictionary to store valid file paths for parameters.
        param (str): The key in `dico_config` corresponding to the file path to validate.
        dico_config (dict): A dictionary containing configuration values, where each key corresponds 
                             to a parameter and its associated value.

    Returns:
        bool: `True` if the file exists and is valid, `False` otherwise.
    """
    link = dico_config[param]
    valid_parameter = string_checking(link)
    if valid_parameter == True:
        if os.path.isfile(link) == False:
            print (BRIGHT_RED + "ERROR... THE FILE " + RESET + link + BRIGHT_RED + " DOES NOT EXIST ..." + RESET)
            return False
        else:
            dico_variables[param] = link
            return True
    else:
        return False


def is_convertible_to_number(s):
    """
    Checks if a string can be converted to a numeric value (integer or float).

    This function attempts to convert the input string `s` to a floating-point number. 
    If the conversion is successful, the function returns `True`, indicating that the 
    string can be interpreted as a number. If the conversion fails (raises a ValueError), 
    the function returns `False`, indicating that the string cannot be converted to a number.

    Args:
        s (str): The string to check for numeric conversion.

    Returns:
        bool: `True` if the string can be converted to a number (int or float), `False` otherwise.
    """
    try:
        float(s)  # try to convert to a float (works for both integers and floats)
        return True
    except ValueError:
        return False


def validate_number(dico_variables : dict, param : str, dico_config : dict) -> bool:
    """
    Validates the numerical value of a specified parameter from the configuration.

    This function checks whether the value of a given parameter in the configuration 
    (`dico_config`) is convertible to a number. It allows the value to be either an 
    integer or a string that can be converted to a number. If the value is valid, it is 
    converted to a float and stored in `dico_variables`. If the value is invalid, an error 
    message is printed, and the function returns `False`.

    Args:
        dico_variables (dict): The dictionary where the validated number will be stored.
        param (str): The name of the parameter to validate.
        dico_config (dict): The dictionary containing the configuration values.

    Returns:
        bool: `True` if the parameter's value is a valid number, `False` if the value is invalid.

    Side Effects:
        - Updates `dico_variables[param]` with the validated number (converted to float).
        - Prints an error message if the value is not a valid number.
    """
    value_param = dico_config[param]
    if type(value_param) == int or type(value_param) == str:
        convert = is_convertible_to_number(value_param)
        if convert == False:
            print (BRIGHT_RED + "ERROR... THE TPM VALUE " + RESET + value_param + BRIGHT_RED + " MUST BE A NUMBER ..." + RESET)
            return False 
        else:
            value_param = float(value_param)
            dico_variables[param] = value_param
            return True
    else:
        print (BRIGHT_RED + "ERROR... THE TPM VALUE " + RESET + str(value_param) + BRIGHT_RED + " MUST BE A NUMBER ..." + RESET)
        return False      
        

def validate_overlap(dico_variables : dict, param : str, dico_config : dict) -> bool:
    """
    Validates the overlap parameter for transcript types.

    This function ensures that the `transcript_overlap` parameter in the configuration 
    is either a string or a list of strings, with each string being one of the accepted values: 
    "intergenic", "intronic", "genic", or "antisense". If the parameter is invalid, an error 
    message is printed. If the parameter is valid, the function updates the `dico_variables[param]` 
    with the validated list of overlap values.

    Args:
        dico_variables (dict): The dictionary where the validated overlap values will be stored.
        param (str): The name of the overlap parameter to validate.
        dico_config (dict): The dictionary containing the configuration values to be validated.

    Returns:
        bool: `True` if the overlap values are valid, `False` if the validation fails.

    Side Effects:
        - Updates `dico_variables[param]` with the validated overlap values.
        - Prints an error message if the overlap value is not valid.
    """
    value_param = dico_config[param]
    if type(value_param) == int or type(value_param) == float:
        print (BRIGHT_RED + "ERROR... THE PARAMETER " + RESET + str(param) + BRIGHT_RED + " MUST NOT BE A NUMBER ..." + RESET)
        return False
    else:
        if type(value_param) == str:
            value_param = [value_param]
        if type(value_param) == list:
            accepted_parameters = ["intergenic", "intronic", "genic", "antisense"]
            correct_param = True
            for i in value_param:
                if i not in accepted_parameters:
                    print (BRIGHT_RED + "ERROR... THE OVERLAP " + RESET + i + BRIGHT_RED + " IS NOT ACCEPTED ..." + RESET)
                    correct_param = False
            if correct_param == False:
                print (BRIGHT_RED + "ACCEPTED OVERLAPS :" + RESET + BRIGHT_GREEN + " intergenic, intronic, genic, antisense" + RESET)
                return False
            else:
                dico_variables[param] = value_param
                return True


def validate_orf_choice(dico_variables : dict, param : str, dico_config : dict) -> bool:
    """
    Validates the ORF (Open Reading Frame) choice parameter.

    This function checks if the specified ORF choice is a valid string. Valid options are 
    "all", "longest", "kozac_highest", and "start_first". If the value is valid, it updates the 
    corresponding entry in the `dico_variables[param][0]` list. If the value is invalid, an error 
    message is printed and the function returns `False`.

    Args:
        dico_variables (dict): The dictionary where the validated ORF choice will be stored.
        param (str): The name of the ORF choice parameter to validate.
        dico_config (dict): The dictionary containing the configuration values to be validated.

    Returns:
        bool: `True` if the ORF choice is a valid string from the allowed options,
              `False` if the validation fails.

    Side Effects:
        - Updates the `dico_variables[param][0]` list with the validated ORF choice value.
        - Prints an error message if the ORF choice is not valid or not a string.
    """
    value_param = dico_config[param]
    if type(value_param) == str:
        if value_param == "all" or value_param == "longest" or value_param == "kozac_highest" or value_param == "start_first":
            dico_variables[param][0] = [value_param]
            return True
        else:
            print (BRIGHT_RED + "ERROR... THE ORF CHOICE " + RESET + value_param + BRIGHT_RED + " IS INVALID ..." + RESET)
            print (BRIGHT_RED + "VALID OPTIONS :" + RESET + BRIGHT_GREEN + " all, longest, kozac_highest, start_first" + RESET)
            return False
    else:
        print (BRIGHT_RED + "ERROR... THE ORF CHOICE MUST BE A STRING" + RESET )
        return False


def validate_utr(dico_variables : dict, param : str, dico_config : dict) -> bool:
    """
    Validates the UTR (Untranslated Region) minimum size parameter.

    This function checks if the specified UTR parameter (either "five_prime" or "three_prime") is an integer
    and whether its value is within the acceptable range (0 to 100, inclusive). If the value is valid, 
    it updates the corresponding value in the `dico_variables["ORFs_choice"]` list. If the value is invalid, 
    an error message is printed and the function returns `False`.

    Args:
        dico_variables (dict): The dictionary where the validated UTR parameter will be stored.
        param (str): The name of the UTR parameter to validate ("five_prime" or "three_prime").
        dico_config (dict): The dictionary containing the configuration values to be validated.

    Returns:
        bool: `True` if the UTR value is a valid integer within the range of 0 to 100, 
              `False` if the validation fails.

    Side Effects:
        - Updates the `dico_variables["ORFs_choice"]` list with the validated UTR minimum size value.
        - Prints an error message if the UTR parameter is not a valid integer or is out of the valid range.
    """
    value_param = dico_config[param]
    if type(value_param) != int:
        print (BRIGHT_RED + "ERROR... THE UTR MIN SIZE " + RESET + str(value_param) + BRIGHT_RED + " MUST BE AN INTEGER ..." + RESET)
        return False
    else:
        if value_param > 100:
            print (BRIGHT_RED + "ERROR... THE UTR MIN SIZE " + RESET + str(value_param) + BRIGHT_RED + " MUST BE BELOW 100 ..." + RESET)
            return False
        else:
            if param == "five_prime":
                dico_variables["ORFs_choice"][2][1] = value_param
            else:
                dico_variables["ORFs_choice"][2][2] = value_param
            return True


def validate_genic_filter(dico_variables : dict, param : str, dico_config : dict) -> bool:
    """
    Validates the genic filter parameter that is represented as a string.

    This function checks if the specified genic filter parameter in the configuration 
    is either the string "True" or "False". If the value is valid, it updates the 
    `dico_variables` dictionary with the corresponding boolean value (`True` or `False`). 
    If the value is invalid, an error message is printed and the function returns `False`.

    Args:
        dico_variables (dict): The dictionary where the validated configuration parameter will be stored.
        param (str): The name of the genic filter parameter to validate.
        dico_config (dict): The dictionary containing the configuration values to be validated.

    Returns:
        bool: `True` if the genic filter value is either "True" or "False" (as strings),
              `False` if the validation fails.

    Side Effects:
        - Updates `dico_variables["filter_genic"]` with the boolean value (`True` or `False`).
        - Prints an error message if the genic filter parameter is not a valid boolean string.
    """
    value_param = dico_config[param]
    if value_param == "True":
        dico_variables["filter_genic"] = True
        return True
    elif value_param == "False":
        dico_variables["filter_genic"] = False
        return True
    else:
        print (BRIGHT_RED + "ERROR... GENIC FILTER " + RESET + str(value_param) + BRIGHT_RED + " MUST BE TRUE OR FALSE ..." + RESET)
        return False


def validate_string_bool(dico_variables : dict, param : str, dico_config : dict) -> bool:
    """
    Validates a boolean parameter that is represented as a string.

    This function checks if the specified parameter in the configuration is either 
    the string "True" or the string "False". If the value is valid, it updates the 
    `dico_variables` dictionary with the validated value. If invalid, an error message 
    is printed and the function returns `False`.

    Args:
        dico_variables (dict): The dictionary where the validated configuration parameter will be stored.
        param (str): The name of the parameter to validate.
        dico_config (dict): The dictionary containing the configuration values to be validated.

    Returns:
        bool: `True` if the parameter is either "True" or "False" (as strings), 
              `False` if the validation fails.
    """
    value_param = dico_config[param]
    if value_param == "True" or value_param == "False":
        dico_variables[param] = value_param
        return True
    else:
        print (BRIGHT_RED + "ERROR... " + str(param) + " " + RESET + " " + str(value_param) + BRIGHT_RED + " MUST BE TRUE OR FALSE ..." + RESET)
        return False


def validate_param_prot_blast(dico_variables : dict, param : str, dico_config : dict) -> bool:
    """
    Validates the 'mode' parameter for a protein BLAST (blastp) configuration.

    This function checks if the 'mode' parameter in the configuration file is one of 
    the acceptable options for protein BLAST. The valid modes are "--sensitive", 
    "--more-sensitive", "--very-sensitive", and "--ultra-sensitive". If the value is 
    valid, it updates the `dico_variables` dictionary with the validated 'mode'. If 
    invalid, an error message is printed and the function returns `False`.

    Args:
        dico_variables (dict): The dictionary holding validated configuration parameters.
        param (str): The parameter to validate (in this case, 'mode').
        dico_config (dict): The dictionary containing the configuration values to be validated.

    Returns:
        bool: `True` if the 'mode' parameter is one of the valid values 
              ("--sensitive", "--more-sensitive", "--very-sensitive", "--ultra-sensitive").
              `False` if the validation fails.
    """
    value_param = dico_config[param]
    if value_param == "--sensitive" or value_param == "--more-sensitive" or value_param == "--very-sensitive" or value_param == "--ultra-sensitive":
        dico_variables[param]["mode"] = value_param
        return True
    else:
        print (BRIGHT_RED + "ERROR... PARAMETERS DIAMOND blastp MUST BE" + RESET )
        print (BRIGHT_GREEN + "--sensitive, --more-sensitive, --very-sensitive or --ultra-sensitive" + RESET )
        return False


def validate_param_nucl_blast(dico_variables : dict, param : str, dico_config : dict) -> bool:
    """
    Validates the 'e_value' parameter for a nucleotide BLAST (blastn) configuration.

    This function checks if the 'e_value' parameter in the configuration file is one 
    of the acceptable values for nucleotide BLAST. The valid options are "0.1", "0.01", 
    "0.001", and "0.0001". If the value is valid, it updates the `dico_variables` dictionary 
    with the validated 'e_value'. If invalid, an error message is printed and the function 
    returns `False`.

    Args:
        dico_variables (dict): The dictionary holding validated configuration parameters.
        param (str): The parameter to validate (in this case, 'e_value').
        dico_config (dict): The dictionary containing the configuration values to be validated.

    Returns:
        bool: `True` if the 'e_value' parameter is one of the valid values ("0.1", "0.01", 
              "0.001", "0.0001"). `False` if the validation fails.
    """
    value_param = dico_config[param]
    if value_param == "0.1" or value_param == "0.01" or value_param == "0.001" or value_param == "0.0001":
        dico_variables[param]["e_value"] = value_param
        return True
    else:
        print (BRIGHT_RED + "ERROR... PARAMETERS DIAMOND blastp MUST BE" + RESET )
        print (BRIGHT_GREEN + "0.1, 0.01, 0.001 or 0.0001" + RESET )
        return False


def validate_synteny_window(dico_variables : dict, param : str, dico_config : dict) -> bool:
    """
    Validates the value of the 'synteny_window' parameter in the configuration.

    This function checks whether the 'synteny_window' parameter from the configuration
    file is an integer within the valid range of 0 to 5, inclusive. If the value is 
    valid, it updates the `dico_variables` dictionary with the validated value. If the 
    value is invalid, it prints an error message and returns `False`.

    Args:
        dico_variables (dict): The dictionary holding validated configuration parameters.
        param (str): The parameter to validate (in this case, "synteny_window").
        dico_config (dict): The dictionary containing the configuration values to be validated.

    Returns:
        bool: `True` if the 'synteny_window' value is a valid integer between 0 and 5, 
              inclusive. `False` if the validation fails (either due to an invalid type or 
              an out-of-range value).
    """
    value_param = dico_config[param]
    if type(value_param) == int:
        if value_param >=0 and value_param <= 5:
            dico_variables[param] = value_param
            return True
        else:
            print (BRIGHT_RED + "ERROR... SYNTENY WINDOW " + RESET + str(value_param) + BRIGHT_RED + " MUST BE BETWEEN 0 AND 5 (both included)" + RESET)
            return False
    else:
        print (BRIGHT_RED + "ERROR... SYNTENY WINDOW " + RESET + str(value_param) + BRIGHT_RED + " MUST BE AN INTEGER" + RESET)
        return False


def validate_premature_stop(dico_variables : dict, param : str, dico_config : dict) -> bool:
    """
    Validates the value of the 'premature_stop' parameter in the configuration.

    This function checks whether the 'premature_stop' parameter from the configuration
    file is an integer within the valid range of 0 to 100, inclusive. If the value is 
    valid, it updates the `dico_variables` dictionary with the validated value. If the 
    value is invalid, it prints an error message and returns `False`.

    Args:
        dico_variables (dict): The dictionary holding validated configuration parameters.
        param (str): The parameter to validate (in this case, "premature_stop").
        dico_config (dict): The dictionary containing the configuration values to be validated.

    Returns:
        bool: `True` if the 'premature_stop' value is a valid integer between 0 and 100, 
              inclusive. `False` if the validation fails (either due to an invalid type or 
              an out-of-range value).
    """
    value_param = dico_config[param]
    if type(value_param) == int:
        if value_param >=0 and value_param <= 100:
            dico_variables[param] = value_param
            return True
        else:
            print (BRIGHT_RED + "ERROR... PREMATURE STOP VALUE " + RESET + str(value_param) + BRIGHT_RED + " MUST BE BETWEEN 0 AND 100 % (both included)" + RESET)
            return False
    else:
        print (BRIGHT_RED + "ERROR... PREMATURE STOP VALUE " + RESET + str(value_param) + BRIGHT_RED + " MUST BE AN INTEGER" + RESET)
        return False


def my_config_file_extract_parameters(link_config : str, strategy : int) -> None:
    """
    Extracts and validates configuration parameters from the provided config file.

    This function reads the configuration file, extracts parameters, and validates 
    each parameter based on the provided strategy (1 or 2). The parameters are 
    stored in a dictionary, `dico_variables`, which is updated with the validated values.
    
    The function checks if required parameters are present in the configuration file 
    and whether their values are of the correct type or format. If any validation fails, 
    the function prints an error message and sets a validation flag to `False`.

    Args:
        link_config (str): The path to the configuration file.
        strategy (int): The strategy to be used, either 1 or 2. Strategy 1 includes 
                        additional parameters (`rec_best_hit`, `synteny_window`, 
                        `premature_stop`), while Strategy 2 omits them.

    Returns:
        tuple: A tuple containing:
            - `dico_variables` (dict): A dictionary of validated configuration parameters.
            - `validate_input` (bool): A flag indicating whether the configuration is valid (`True`) or not (`False`).
        
    Side Effects:
        - Updates `dico_variables` with validated values.
        - Prints error messages if any validation fails, indicating which parameter is invalid or missing.
    """
    dico_config = parse_config(link_config)
    dico_variables = {"strategy": "1",
                        "query" : "",
                        "path_to_genome_repository":"", 
                        "path_to_transcriptome_repository":"", 
                        "link_database_outgroup_prot":"", 
                        "link_database_outgroup_nucl":"", 
                        "TPM_threeshold": 0.5,
                        "transcript_overlap" : ["intergenic"],
                        "ORFs_choice" : [["longest"],["duplicate_handle"], ["utr_size",0,0]],
                        "filter_genic" : False,
                        "filter_TE" : "False",
                        "rm_undir_transc" : "False",
                        "parameters_database_prot" : {"type" : "blastp", "mode" : "--more-sensitive"},
                        "parameters_database_nucl" : {"type" : "blastn", "e_value" : "0.01", "coverage" : None, "strand" : None},
                        "rec_best_hit" : "False",
                        "synteny_window" : 2,
                        "premature_stop": 50}

    if strategy == 2:
        dico_variables["strategy"] = 2
        del dico_variables["rec_best_hit"]
        del dico_variables["synteny_window"]
        del dico_variables["premature_stop"]

    validate_input = True
    if "query" in dico_config:
        if type(dico_config["query"]) == str:
            dico_variables["query"] = dico_config["query"]
            validate_input = True
        else:
            print (BRIGHT_RED + "ERROR... THE QUERY NAME " + RESET + str(dico_config["query"]) + BRIGHT_RED + " MUST BE A STRING ..." + RESET)
            validate_input = False
    else:
        print (BRIGHT_RED + "ERROR... THE QUERY NAME IS MISSING" + RESET)
        validate_input = False
    if validate_input == True and "path_to_genome_repository" in dico_config:
        validate_input = check_true_repository(dico_variables, "path_to_genome_repository", dico_config)
    if validate_input == True and "path_to_transcriptome_repository" in dico_config:
        validate_input = check_true_repository(dico_variables, "path_to_transcriptome_repository", dico_config)
    if validate_input == True and "link_database_outgroup_prot" in dico_config:
        validate_input = check_true_file(dico_variables, "link_database_outgroup_prot", dico_config)
    if validate_input == True and "link_database_outgroup_nucl" in dico_config:
        validate_input = check_true_file(dico_variables, "link_database_outgroup_nucl", dico_config)
    if validate_input == True and "TPM_threeshold" in dico_config:
        validate_input = validate_number(dico_variables, "TPM_threeshold", dico_config)
    if validate_input == True and "transcript_overlap" in dico_config:
        validate_input = validate_overlap(dico_variables, "transcript_overlap", dico_config)
    if validate_input == True and "ORFs_choice" in dico_config:
        validate_input = validate_orf_choice(dico_variables, "ORFs_choice", dico_config)
    if validate_input == True and "five_prime" in dico_config:
        validate_input = validate_utr(dico_variables, "five_prime", dico_config)
    if validate_input == True and "three_prime" in dico_config:
        validate_input = validate_utr(dico_variables, "three_prime", dico_config)
    if validate_input == True and "filter_genic" in dico_config:
        validate_input = validate_genic_filter(dico_variables, "filter_genic", dico_config)
    if validate_input == True and "filter_TE" in dico_config:
        validate_input = validate_string_bool(dico_variables, "filter_TE", dico_config)
    if validate_input == True and "rm_undir_transc" in dico_config:
        validate_input = validate_string_bool(dico_variables, "rm_undir_transc", dico_config)
    if validate_input == True and "parameters_database_prot" in dico_config:
        validate_input = validate_param_prot_blast(dico_variables, "parameters_database_prot", dico_config)
    if validate_input == True and "parameters_database_nucl" in dico_config:
        validate_input = validate_param_nucl_blast(dico_variables, "parameters_database_nucl", dico_config)
    if strategy == 1:
        if validate_input == True and "rec_best_hit" in dico_config:
            validate_input = validate_string_bool(dico_variables, "rec_best_hit", dico_config)
        if validate_input == True and "synteny_window" in dico_config:
            validate_input = validate_synteny_window(dico_variables, "synteny_window", dico_config)
        if validate_input == True and "premature_stop" in dico_config:
            validate_input = validate_premature_stop(dico_variables, "premature_stop", dico_config)
    return dico_variables, validate_input

