import os
import sys
from deswoman.Main_run_part1 import run_part_1_strat1
from deswoman.Main_run_part1 import run_part_1_strat2
from deswoman.Main_run_part2 import run_part_2
from deswoman.Main_run_part3 import run_part3_strat1
from deswoman.Main_run_part3 import run_part3_strat2
from deswoman.module_colors import *
from deswoman.module_input_checking_strat1 import assess_parameters_strat1
from deswoman.module_input_checking_strat2 import assess_parameters_strat2


def main(strategy: str, link_config: bool | str) -> None:
    """
    Executes the main workflow of DESwoMAN based on the selected strategy.

    This function runs three successive parts of the DESwoMAN process, depending on whether
    "Strategy1" or "Strategy2" is chosen. It first validates user-defined parameters through
    a graphical interface and ensures input correctness before proceeding.

    Args:
        strategy (str): The selected strategy, either "Strategy1" or "Strategy2".
        link_config (bool|str): The link to the config file, or FALSE which indicates DESwoMAN runs with the graphical interface

    Workflow:
        1. Validates parameters through `assess_parameters_strat1` or `assess_parameters_strat2`.
        2. Runs Part 1, which detects ORFs.
        3. If ORFs are found, Part 2 is executed to test homology against an outgroup dataset.
        4. Part 3 is executed based on the selected strategy.

    Notes:
        - If parameter validation fails, execution stops with appropriate error messages.
        - Part 2 is skipped if no outgroup database is provided.
        - The process concludes with a completion message if all parts execute successfully.

    Returns:
        None
    """
    # If the user asks for Strategy 1 or 2, DESwoMAN opens the graphical interface of strat1 or 2, so that the user choses parameters. DESwoMAN then assesses
    # a couple of things to make sure that the input are valids.
    if strategy == "Strategy1":
        valid_parameters, dico_variables = assess_parameters_strat1(link_config)
    else:
        valid_parameters, dico_variables = assess_parameters_strat2(link_config)
    # DESwoMAN only run if all parameters were validated; If not, the previous function display the error messages corresponding to te problem.
    if valid_parameters == True:
        validated_step1 = False
        print(BRIGHT_GREEN + "- - - - - - - - - - - - - - - - - - -" + RESET)
        print("")
        print(BRIGHT_YELLOW + "Running DESwoMAN Part 1." + RESET)

        if dico_variables["strategy"] == "1":
            # Runs Part 1 Strategy 1
            validated_step1 = run_part_1_strat1(dico_variables)
        else:
            # Runs Part 1 Strategy 2
            validated_step1 = run_part_1_strat2(dico_variables)

        print(YELLOW + "DESwoMAN Part 1 : DONE!" + RESET)
        print("")
        # if no ORFs were detected, then DESwoMAN stops after part 1 as the variable "validated_step1" becomes False
        if validated_step1 == True:
            # Part 2 runs only if the user specified at least one outgroup dataset to test homology
            if (
                dico_variables["link_database_outgroup_prot"] != ""
                or dico_variables["link_database_outgroup_nucl"] != ""
            ):
                print(BRIGHT_YELLOW + "Running DESwoMAN Part 2." + RESET)
                print("")
                run_part_2(dico_variables)
                print(YELLOW + "DESwoMAN Part 2 : DONE!" + RESET)
            else:
                print(
                    BRIGHT_YELLOW
                    + "DESwoMAN Part 2 NOT PERFORMED (no database for homology search)"
                    + RESET
                )
            # Start part 3
            if dico_variables["strategy"] == "1":
                print("")
                print(BRIGHT_YELLOW + "Running DESwoMAN Part 3." + RESET)
                run_part3_strat1(dico_variables)
            else:
                print("")
                print(BRIGHT_YELLOW + "Running DESwoMAN Part 3." + RESET)
                run_part3_strat2(dico_variables)

            # End of the RUN
            print(YELLOW + "DESwoMAN Part 3 : DONE!" + RESET)
            print(" ")
            print(BRIGHT_YELLOW + "Goodbye ;) !" + RESET)


def validate_link_config(link: str) -> bool:
    """
    Validates the existence of a configuration file at the given file path.

    This function checks if the file specified by the `link` exists on the file system.
    If the file exists, it returns `True`. If the file does not exist, it prints an error
    message and returns `False`.

    Args:
        link (str): The path to the configuration file to be validated.

    Returns:
        bool: `True` if the file exists, `False` otherwise.

    Side Effects:
        - Prints an error message if the file does not exist, including the incorrect file path.
    """
    if os.path.isfile(link) == False:
        print(
            BRIGHT_RED
            + "ERROR : The link to config file is not correct : "
            + RESET
            + link
        )
        return False
    else:
        return True


def entry_point_main():
    """
    Main entry point for running the DESwoMAN script.

    This block handles command-line arguments to determine the strategy to run
    and to provide the necessary configuration file link. It also performs error
    handling for incorrect or missing arguments.

    Command-line Arguments:
        - strategy (str): The first argument specifies the strategy to use,
                          either "Strategy1" or "Strategy2".
        - link_config (str, optional): The second argument (if provided) is the
                                       path to the configuration file. It is
                                       validated before being passed to the main function.

    Workflow:
        1. If a strategy argument is provided, it checks if it's "Strategy1" or "Strategy2".
        2. If a valid strategy is provided and a configuration file is given, it validates the link
           to the configuration file.
        3. The `main()` function is called with the strategy and configuration, if applicable.
        4. If the arguments are invalid (wrong strategy or missing configuration),
           appropriate error messages are displayed.

    Error Handling:
        - If no strategy is provided or an invalid strategy is entered, an error message is displayed.
        - If the configuration file link is invalid or missing, an error message is displayed.

    Returns:
        None
    """
    if len(sys.argv) > 1:
        strategy = sys.argv[1]  # The first argument after the script name
        if strategy == "Strategy1" or strategy == "Strategy2":
            CONFIG = False
            if len(sys.argv) > 2:
                link_config = sys.argv[2]
                CONFIG = validate_link_config(link_config)
                if CONFIG == True:
                    main(strategy, link_config)
            else:
                main(strategy, CONFIG)
        else:
            print(
                BRIGHT_RED
                + "ERROR : The chosen strategy is not correct : pick Strategy1 or Strategy2"
                + RESET
            )
    else:
        print(
            BRIGHT_RED
            + "ERROR : No strategy was entered. A strategy must be entered in order to run DESwoMAN"
            + RESET
        )


if __name__ == "__main__":
    entry_point_main()
