from tkinter import *
import customtkinter
from tkinter.filedialog import *
from tkinter import filedialog


__author__ = "Anna Grandchamp"
__contributor__ = ""
__copyright__ = ""
__credits__ = []
__license__ = ""
__version__ = "1.0.0"
__maintainer__ = "Anna Grandchamp"
__email__ = "anna.grandchamp@inserm.fr"


###########################
# window top left : Input data
###########################


def graphical_select_genome_folder_strat1(
    dico_file_and_dir: dict,
    radiobutton_frame: customtkinter.CTkFrame,
    win: customtkinter.CTk,
    value_ref_name: StringVar,
) -> None:
    """
    Opens a file dialog to select a genome folder and updates the GUI accordingly.

    This function allows the user to select a directory for genome data storage.
    It updates the `dico_file_and_dir` dictionary with the selected path and
    displays a confirmation label in the provided GUI frame. If both the
    transcriptome repository path and reference name are set, it enables a "RUN" button.

    Args:
        dico_file_and_dir (dict): A dictionary storing file paths and directory selections.
            - Updates the "path_to_genome_repository" key with the selected folder path.
        radiobutton_frame (customtkinter.CTkFrame): The GUI frame where the confirmation label is displayed.
        win (customtkinter.CTk): The main application window, used for placing the "RUN" button.
        value_ref_name (tk.StringVar): A Tkinter string variable containing the reference name.

    Side Effects:
        - Opens a folder selection dialog for the user.
        - Updates `dico_file_and_dir` with the selected genome folder path.
        - Displays a confirmation label ("V") if a folder is selected.
        - Enables a "RUN" button if both the transcriptome repository and reference name are set.

    Dependencies:
        - Uses `filedialog.askdirectory()` to prompt folder selection.
        - Relies on `customtkinter` for UI elements.

    """
    folder_selected = filedialog.askdirectory()
    dico_file_and_dir["path_to_genome_repository"] = folder_selected
    label = customtkinter.CTkLabel(
        radiobutton_frame,
        text="V",
        width=30,
        height=25,
        text_color="#37589c",
        corner_radius=4,
        fg_color="white",
        font=("Arial", 15),
    ).place(relx=0.38, rely=0.61)
    name_ref = value_ref_name.get()
    if dico_file_and_dir["path_to_transcriptome_repository"] != "" and name_ref != "":
        button_DESMAN_run = customtkinter.CTkButton(
            master=win,
            width=100,
            height=50,
            text="RUN",
            command=lambda: graphical_Close_strat1(win),
        ).place(relx=0.05, rely=0.91, anchor="w")


def graphical_select_transcriptome_folder_strat1(
    dico_file_and_dir: dict,
    radiobutton_frame: customtkinter.CTkFrame,
    win: Tk,
    value_ref_name: StringVar,
) -> None:
    """
    Opens a file dialog to select a transcriptome folder and updates the GUI accordingly.

    This function allows the user to select a directory for transcriptome data storage.
    It updates the `dico_file_and_dir` dictionary with the selected path and
    displays a confirmation label in the provided GUI frame. If both the
    genome repository path and reference name are set, it enables a "RUN" button.

    Args:
        dico_file_and_dir (dict): A dictionary storing file paths and directory selections.
            - Updates the "path_to_transcriptome_repository" key with the selected folder path.
        radiobutton_frame (customtkinter.CTkFrame): The GUI frame where the confirmation label is displayed.
        win (Tk): The main application window, used for placing the "RUN" button.
        value_ref_name (StringVar): A Tkinter string variable containing the reference name.

    Side Effects:
        - Opens a folder selection dialog for the user.
        - Updates `dico_file_and_dir` with the selected transcriptome folder path.
        - Displays a confirmation label ("V") if a folder is selected.
        - Enables a "RUN" button if both the genome repository and reference name are set.

    Dependencies:
        - Uses `filedialog.askdirectory()` to prompt folder selection.
        - Relies on `customtkinter` for UI elements.

    Returns:
        None
    """
    folder_selected = filedialog.askdirectory()
    dico_file_and_dir["path_to_transcriptome_repository"] = folder_selected
    label = customtkinter.CTkLabel(
        radiobutton_frame,
        text="V",
        width=30,
        height=25,
        text_color="#37589c",
        corner_radius=4,
        fg_color="white",
        font=("Arial", 15),
    ).place(relx=0.905, rely=0.61)
    name_ref = value_ref_name.get()
    if dico_file_and_dir["path_to_genome_repository"] != "" and name_ref != "":
        button_DESMAN_run = customtkinter.CTkButton(
            master=win,
            width=100,
            height=50,
            text="RUN",
            command=lambda: graphical_Close_strat1(win),
        ).place(relx=0.05, rely=0.91, anchor="w")


def graphical_Close_strat1(win: Tk) -> None:
    """
    Closes the given Tkinter window.

    This function destroys the specified Tkinter window, effectively closing the GUI.

    Args:
        win (Tk): The Tkinter window instance to be closed.

    Returns:
        None
    """
    win.destroy()


def graphical_activate_with_text_strat1(win: Tk, dico_file_and_dir: dict) -> None:
    """
    Activates the "RUN" button if both genome and transcriptome repository paths are set.

    This function checks if both the "path_to_genome_repository" and "path_to_transcriptome_repository"
    keys in the `dico_file_and_dir` dictionary are not empty. If so, it creates and displays a "RUN"
    button on the provided Tkinter window, which will trigger the closure of the window when clicked.

    Args:
        win (Tk): The Tkinter window where the "RUN" button will be placed.
        dico_file_and_dir (dict): A dictionary storing paths to the genome and transcriptome repositories.
            - Keys "path_to_genome_repository" and "path_to_transcriptome_repository" should be non-empty
              to enable the "RUN" button.

    Returns:
        None
    """
    if (
        dico_file_and_dir["path_to_genome_repository"] != ""
        and dico_file_and_dir["path_to_transcriptome_repository"] != ""
    ):
        button_DESMAN_run = customtkinter.CTkButton(
            master=win,
            width=100,
            height=50,
            text="RUN",
            command=lambda: graphical_Close_strat1(win),
        ).place(relx=0.05, rely=0.91, anchor="w")


def graphical_enter_mandatory_strat1(
    win: Tk,
    radiobutton_frame: customtkinter.CTkFrame,
    my_font: tuple,
    dico_file_and_dir: dict,
) -> StringVar:
    """
    Displays the mandatory fields for selecting genome and transcriptome folders and entering a query name.

    This function sets up the necessary graphical elements for entering a query name,
    selecting the genome directory, and selecting the transcriptome directory. It binds
    the query name entry field to enable the "RUN" button once both the genome and
    transcriptome repositories are selected. It returns the StringVar associated with the query name.

    Args:
        win (Tk): The Tkinter window where the graphical elements will be placed.
        radiobutton_frame (customtkinter.CTkFrame): The frame in which the labels and buttons are placed.
        my_font (tuple): The font configuration for the labels.
        dico_file_and_dir (Dict[str, str]): A dictionary storing paths to the genome and transcriptome directories.

    Returns:
        StringVar: The StringVar instance containing the query name entered by the user.
    """
    # query name
    label = customtkinter.CTkLabel(
        radiobutton_frame, text="Enter query name : ", font=my_font
    ).place(relx=0.02, rely=0.3, anchor="w")
    value_ref_name = StringVar(value="")
    entree = Entry(radiobutton_frame, textvariable=value_ref_name)
    entree.place(relx=0.40, rely=0.28, anchor="w")
    entree.bind(
        "<KeyRelease>",
        lambda e: graphical_activate_with_text_strat1(win, dico_file_and_dir),
    )
    # buttons
    label = customtkinter.CTkLabel(
        radiobutton_frame,
        text="",
        font=my_font,
        corner_radius=4,
        width=30,
        height=25,
        fg_color="white",
    ).place(relx=0.38, rely=0.61)
    label = customtkinter.CTkLabel(
        radiobutton_frame,
        text="",
        font=my_font,
        corner_radius=4,
        width=30,
        height=25,
        fg_color="white",
    ).place(relx=0.905, rely=0.61)
    button_genome_dir = customtkinter.CTkButton(
        master=radiobutton_frame,
        width=150,
        height=25,
        text="Genome directory",
        command=lambda: graphical_select_genome_folder_strat1(
            dico_file_and_dir, radiobutton_frame, win, value_ref_name
        ),
    ).place(relx=0.02, rely=0.75, anchor="w")
    button_transcriptome_dir = customtkinter.CTkButton(
        master=radiobutton_frame,
        width=150,
        height=25,
        text="Transcriptome directory",
        command=lambda: graphical_select_transcriptome_folder_strat1(
            dico_file_and_dir, radiobutton_frame, win, value_ref_name
        ),
    ).place(relx=0.54, rely=0.75, anchor="w")
    return value_ref_name


###########################
### dataset for Step 2 Blast homology search
###########################


def graphical_select_prot_database_strat1(
    sidebar_frame_left: customtkinter.CTkFrame, dico_file_and_dir: dict
) -> None:
    """
    Opens a file dialog to select a protein database for homology search and updates the GUI accordingly.

    This function allows the user to select a protein dataset file (e.g., a FASTA file) and stores the
    selected file path in the provided `dico_file_and_dir` dictionary under the key "link_database_outgroup_prot".
    It also displays a confirmation label in the provided sidebar frame to indicate that a dataset has been selected.

    Args:
        sidebar_frame_left (customtkinter.CTkFrame): The GUI frame where the confirmation label is displayed.
        dico_file_and_dir (Dict[str, str]): A dictionary storing file paths and other directory-related information.
            - The key "link_database_outgroup_prot" will be updated with the selected file path.

    Returns:
        None
    """
    filepath = askopenfilename(
        title="Select a dataset",
        filetypes=[("fasta file", ".fa .fna .fasta"), ("all files", ".*")],
    )
    dico_file_and_dir["link_database_outgroup_prot"] = filepath
    label = customtkinter.CTkLabel(
        sidebar_frame_left,
        text="V",
        width=30,
        height=25,
        text_color="#37589c",
        corner_radius=4,
        font=("Arial", 15),
    ).place(relx=0.36, rely=0.418)


def graphical_select_nucl_database_strat1(
    sidebar_frame_left: customtkinter.CTkFrame, dico_file_and_dir: dict
) -> None:
    """
    Opens a file dialog to select a nucleotide database for homology search and updates the GUI accordingly.

    This function allows the user to select a nucleotide dataset file (e.g., a FASTA file) and stores the
    selected file path in the provided `dico_file_and_dir` dictionary under the key "link_database_outgroup_nucl".
    It also displays a confirmation label in the provided sidebar frame to indicate that a dataset has been selected.

    Args:
        sidebar_frame_left (customtkinter.CTkFrame): The GUI frame where the confirmation label is displayed.
        dico_file_and_dir (Dict[str, str]): A dictionary storing file paths and other directory-related information.
            - The key "link_database_outgroup_nucl" will be updated with the selected file path.

    Returns:
        None
    """
    filepath = askopenfilename(
        title="Select a dataset",
        filetypes=[("fasta file", ".fa .fna .fasta"), ("all files", ".*")],
    )
    dico_file_and_dir["link_database_outgroup_nucl"] = filepath
    label = customtkinter.CTkLabel(
        sidebar_frame_left,
        text="V",
        width=30,
        height=25,
        text_color="#37589c",
        corner_radius=4,
        font=("Arial", 15),
    ).place(relx=0.36, rely=0.72)


def graphical_enter_databases_strat1(
    sidebar_frame_left: customtkinter.CTkFrame, my_font: tuple, dico_file_and_dir: dict
) -> tuple:
    """
    Handles the selection of BLAST datasets and sets up the options for protein and nucleotide databases.

    This function creates buttons to select a protein or nucleotide BLAST database,
    and it allows the user to select the sensitivity mode for BLAST searches. The function
    also returns the selected sensitivity options for protein and nucleotide searches.

    Args:
        sidebar_frame_left (customtkinter.CTkFrame): The GUI frame where the buttons and radio buttons are placed.
        my_font (tuple): The font configuration for the labels and buttons.
        dico_file_and_dir (Dict[str, str]): A dictionary storing file paths and other directory-related information.
            - The keys "link_database_outgroup_prot" and "link_database_outgroup_nucl" will be updated.

    Returns:
        Tuple[IntVar, IntVar]: Returns two `IntVar` instances, one for the protein BLAST sensitivity mode
        and one for the nucleotide BLAST sensitivity mode.
    """
    button_prot_file = customtkinter.CTkButton(
        master=sidebar_frame_left,
        width=150,
        height=25,
        text="BLAST protein database",
        command=lambda: graphical_select_prot_database_strat1(
            sidebar_frame_left, dico_file_and_dir
        ),
    ).place(relx=0.02, rely=0.46, anchor="w")
    button_nucl_file = customtkinter.CTkButton(
        master=sidebar_frame_left,
        width=150,
        height=25,
        text="BLAST nucl database",
        command=lambda: graphical_select_nucl_database_strat1(
            sidebar_frame_left, dico_file_and_dir
        ),
    ).place(relx=0.02, rely=0.76, anchor="w")
    var_blastp_mode = IntVar(value=2)
    # Define a Checkbox
    case_1 = customtkinter.CTkRadioButton(
        master=sidebar_frame_left, text="sensitive", variable=var_blastp_mode, value=1
    ).place(relx=0.02, rely=0.58, anchor="w")
    case_2 = customtkinter.CTkRadioButton(
        master=sidebar_frame_left,
        text="more sensitive",
        variable=var_blastp_mode,
        value=2,
    ).place(relx=0.22, rely=0.58, anchor="w")
    case_3 = customtkinter.CTkRadioButton(
        master=sidebar_frame_left,
        text="very sensitive",
        variable=var_blastp_mode,
        value=3,
    ).place(relx=0.49, rely=0.58, anchor="w")
    case_4 = customtkinter.CTkRadioButton(
        master=sidebar_frame_left,
        text="ultra sensitive",
        variable=var_blastp_mode,
        value=4,
    ).place(relx=0.75, rely=0.58, anchor="w")
    var_blastnucl_mode = IntVar(value=2)
    # Define a Checkbox
    case_5 = customtkinter.CTkRadioButton(
        master=sidebar_frame_left, text="0.1", variable=var_blastnucl_mode, value=1
    ).place(relx=0.02, rely=0.88, anchor="w")
    case_6 = customtkinter.CTkRadioButton(
        master=sidebar_frame_left, text="0.01", variable=var_blastnucl_mode, value=2
    ).place(relx=0.16, rely=0.88, anchor="w")
    case_7 = customtkinter.CTkRadioButton(
        master=sidebar_frame_left, text="0.001", variable=var_blastnucl_mode, value=3
    ).place(relx=0.32, rely=0.88, anchor="w")
    case_8 = customtkinter.CTkRadioButton(
        master=sidebar_frame_left, text="0.0001", variable=var_blastnucl_mode, value=4
    ).place(relx=0.50, rely=0.88, anchor="w")
    return var_blastp_mode, var_blastnucl_mode


###########################
# window rigth : Step 1 optional parameters
###########################


def graphical_TPM_threeshold_strat1(
    sidebar_frame_right: customtkinter.CTkFrame, my_font: tuple
) -> StringVar:
    """
    Allows the user to set the minimum TPM (Transcripts Per Million) threshold for transcripts.

    This function creates a label and a Spinbox widget in the provided sidebar frame that allows the user to
    select the minimum TPM value for transcripts. The value is stored in a `StringVar` instance.

    Args:
        sidebar_frame_right (customtkinter.CTkFrame): The GUI frame where the label and Spinbox are displayed.
        my_font (tuple): The font configuration for the label text.

    Returns:
        StringVar: A `StringVar` instance holding the value selected by the user for the minimum TPM threshold.
    """
    label = customtkinter.CTkLabel(
        sidebar_frame_right, text="Minimum TPM level for transcripts", font=my_font
    ).place(relx=0.02, rely=0.05, anchor="w")
    var_TPM = StringVar(value=0)
    choice_TPM = Spinbox(
        sidebar_frame_right,
        from_=0.5,
        to=10,
        increment=0.5,
        textvariable=var_TPM,
        wrap=True,
    ).place(relx=0.02, rely=0.1, anchor="w")
    return var_TPM


def graphical_filter_TE_strat1(
    sidebar_frame_right: customtkinter.CTkFrame, my_font: tuple
) -> IntVar:
    """
    Allows the user to decide whether to filter TEs (Transposable Elements).

    This function creates a label and two radio buttons (True/False) in the provided sidebar frame,
    enabling the user to select whether or not to filter transposable elements (TEs). The selection
    is stored in an `IntVar` instance.

    Args:
        sidebar_frame_right (customtkinter.CTkFrame): The GUI frame where the label and radio buttons are displayed.
        my_font (tuple): The font configuration for the label text.

    Returns:
        IntVar: An `IntVar` instance holding the value selected by the user for filtering TEs.
    """
    try:
        # Define empty variables
        var_filter_TE = customtkinter.IntVar(value=1)

        label = customtkinter.CTkLabel(
            sidebar_frame_right, text="Filter TEs", font=my_font
        ).place(relx=0.02, rely=0.18, anchor="w")
        # Define a Checkbox
        case_filter_TE = customtkinter.CTkRadioButton(
            master=sidebar_frame_right, text="False", variable=var_filter_TE, value=1
        ).place(relx=0.02, rely=0.22, anchor="w")
        case_no_filter_TEs = customtkinter.CTkRadioButton(
            master=sidebar_frame_right, text="True", variable=var_filter_TE, value=2
        ).place(relx=0.18, rely=0.22, anchor="w")
    except BaseException as e:
        pass
    return var_filter_TE


def graphical_remove_undir_transcripts_strat1(
    sidebar_frame_right: customtkinter.CTkFrame, my_font: tuple
) -> IntVar:
    """
    Allows the user to decide whether to remove transcripts with a "." direction.

    This function creates a label and two radio buttons (True/False) in the provided sidebar frame,
    enabling the user to select whether or not to remove transcripts with the "." direction. The selection
    is stored in an `IntVar` instance.

    Args:
        sidebar_frame_right (customtkinter.CTkFrame): The GUI frame where the label and radio buttons are displayed.
        my_font (tuple): The font configuration for the label text.

    Returns:
        IntVar: An `IntVar` instance holding the value selected by the user for removing unoriented transcripts.
    """
    # Define empty variables
    var_undir_transcr = IntVar(value=1)
    label = customtkinter.CTkLabel(
        sidebar_frame_right, text="Remove unoriented transcripts", font=my_font
    ).place(relx=0.02, rely=0.28, anchor="w")
    # Define a Checkbox
    case_remove_undir_tr = customtkinter.CTkRadioButton(
        master=sidebar_frame_right, text="False", variable=var_undir_transcr, value=1
    ).place(relx=0.02, rely=0.32, anchor="w")
    case_no_remove_undir_tr = customtkinter.CTkRadioButton(
        master=sidebar_frame_right, text="True", variable=var_undir_transcr, value=2
    ).place(relx=0.18, rely=0.32, anchor="w")
    return var_undir_transcr


def graphical_transcript_location_strat1(
    sidebar_frame_right: customtkinter.CTkFrame, my_font: tuple
) -> tuple:
    """
    Allows the user to select which transcripts to keep according to their genomic location.

    This function creates a label and four checkboxes in the provided sidebar frame,
    enabling the user to choose whether to keep transcripts based on their location (intergenic,
    intronic, antisense, or genic). The selections are stored in `IntVar` instances.

    Args:
        sidebar_frame_right (customtkinter.CTkFrame): The GUI frame where the checkboxes are displayed.
        my_font (tuple): The font configuration for the label text.

    Returns:
        Tuple[IntVar, IntVar, IntVar, IntVar]: A tuple containing four `IntVar` instances corresponding
        to the user's selections for each transcript location type (intergenic, intronic, antisense, and genic).
    """
    # Define empty variables
    var_intergenic = IntVar(value=1)
    var_intronic = IntVar()
    var_antisense = IntVar()
    var_genic = IntVar()
    label = customtkinter.CTkLabel(
        sidebar_frame_right, text="Transcript location", font=my_font
    ).place(relx=0.02, rely=0.38, anchor="w")
    # Define a Checkbox
    case_intergenic = customtkinter.CTkCheckBox(
        master=sidebar_frame_right,
        text="intergenic",
        variable=var_intergenic,
        onvalue=1,
        offvalue=0,
    ).place(relx=0.02, rely=0.42, anchor="w")
    case_intronic = customtkinter.CTkCheckBox(
        master=sidebar_frame_right,
        text="intronic",
        variable=var_intronic,
        onvalue=1,
        offvalue=0,
    ).place(relx=0.29, rely=0.42, anchor="w")
    case_antisense = customtkinter.CTkCheckBox(
        master=sidebar_frame_right,
        text="antisense",
        variable=var_antisense,
        onvalue=1,
        offvalue=0,
    ).place(relx=0.56, rely=0.42, anchor="w")
    case_genic = customtkinter.CTkCheckBox(
        master=sidebar_frame_right,
        text="genic",
        variable=var_genic,
        onvalue=1,
        offvalue=0,
    ).place(relx=0.83, rely=0.42, anchor="w")
    return var_intergenic, var_intronic, var_antisense, var_genic


def graphical_ORFs_choice_strat1(
    sidebar_frame_right: customtkinter.CTkFrame, my_font: tuple
) -> IntVar:
    """
    Allows the user to select which ORFs (Open Reading Frames) to consider in transcripts.

    This function creates a label and several radio buttons that enable the user to choose
    between different criteria for selecting ORFs in the transcripts (e.g., all ORFs, the longest ORF,
    the ORF with the highest Kozak score, or the first ORF in the transcript). The user's selection
    is stored in an `IntVar` instance.

    Args:
        sidebar_frame_right (customtkinter.CTkFrame): The GUI frame where the radio buttons are displayed.
        my_font (tuple): The font configuration for the label text.

    Returns:
        IntVar: A variable holding the user's selection for the ORF choice.
               Possible values:
               - 1: All ORFs
               - 2: Longest ORF
               - 3: Highest Kozak score
               - 4: First ORF in transcript
    """
    # Define empty variables
    var_orf_choice = IntVar(value=2)
    label = customtkinter.CTkLabel(
        sidebar_frame_right, text="ORF choice", font=my_font
    ).place(relx=0.02, rely=0.49, anchor="w")
    # Define a Checkbox
    case_all_orfs = customtkinter.CTkRadioButton(
        master=sidebar_frame_right, text="All ORFS", variable=var_orf_choice, value=1
    ).place(relx=0.02, rely=0.53, anchor="w")
    case_intronic = customtkinter.CTkRadioButton(
        master=sidebar_frame_right, text="Longest ORF", variable=var_orf_choice, value=2
    ).place(relx=0.02, rely=0.58, anchor="w")
    case_antisense = customtkinter.CTkRadioButton(
        master=sidebar_frame_right,
        text="Highest Kozac score",
        variable=var_orf_choice,
        value=3,
    ).place(relx=0.02, rely=0.63, anchor="w")
    case_genic = customtkinter.CTkRadioButton(
        master=sidebar_frame_right,
        text="First ORF in transcript",
        variable=var_orf_choice,
        value=4,
    ).place(relx=0.02, rely=0.68, anchor="w")
    return var_orf_choice


def graphical_min_size_5UTR_strat1(
    sidebar_frame_right: customtkinter.CTkFrame, my_font: tuple
) -> StringVar:
    """
    Allows the user to specify the minimum length for the 5'UTR (Untranslated Region) in nucleotides.

    This function creates a label and a Spinbox widget for the user to enter a value for the minimum size
    of the 5'UTR. The user can input values between 0 and 100 nucleotides. The user's input is stored in
    a `StringVar` instance.

    Args:
        sidebar_frame_right (customtkinter.CTkFrame): The GUI frame where the Spinbox widget is displayed.
        my_font (tuple): The font configuration for the label text.

    Returns:
        StringVar: A variable holding the user's input for the minimum 5'UTR size.
    """
    label = customtkinter.CTkLabel(
        sidebar_frame_right, text="Min 5'UTR size (nucl)", font=my_font
    ).place(relx=0.02, rely=0.75, anchor="w")
    var_5UTR = StringVar(value=0)
    choice_5prime_min_size = Spinbox(
        sidebar_frame_right, from_=0, to=100, textvariable=var_5UTR, wrap=True
    ).place(relx=0.02, rely=0.79, anchor="w")
    return var_5UTR


def graphical_min_size_3UTR_strat1(
    sidebar_frame_right: customtkinter.CTkFrame, my_font: tuple
) -> StringVar:
    """
    Allows the user to specify the minimum length for the 3'UTR (Untranslated Region) in nucleotides.

    This function creates a label and a Spinbox widget for the user to enter a value for the minimum size
    of the 3'UTR. The user can input values between 0 and 100 nucleotides. The user's input is stored in
    a `StringVar` instance.

    Args:
        sidebar_frame_right (customtkinter.CTkFrame): The GUI frame where the Spinbox widget is displayed.
        my_font (tuple): The font configuration for the label text.

    Returns:
        StringVar: A variable holding the user's input for the minimum 3'UTR size.
    """
    label = customtkinter.CTkLabel(
        sidebar_frame_right, text="Min 3'UTR size (nucl)", font=my_font
    ).place(relx=0.54, rely=0.75, anchor="w")
    var_3UTR = StringVar(value=0)
    choice_3prime_min_size = Spinbox(
        sidebar_frame_right, from_=0, to=100, textvariable=var_3UTR, wrap=True
    ).place(relx=0.56, rely=0.79, anchor="w")
    return var_3UTR


def graphical_filter_genic_transcripts_strat1(
    sidebar_frame_right: customtkinter.CTkFrame, my_font: tuple
) -> IntVar:
    """
    Allows the user to decide whether to filter splice variants of genic transcripts that correspond to new
    splicing events of known genes, even though they do not overlap with the gene.

    This function creates radio buttons that allow the user to choose whether to filter these transcripts
    (True or False). The user's choice is stored in an `IntVar` instance.

    Args:
        sidebar_frame_right (customtkinter.CTkFrame): The GUI frame where the radio buttons are displayed.
        my_font (tuple): The font configuration for the label text.

    Returns:
        IntVar: A variable holding the user's choice to filter genic transcripts.
    """
    # Define empty variables
    var_filter_genic = IntVar(value=1)
    label = customtkinter.CTkLabel(
        sidebar_frame_right,
        text="Filter splice variants of genic transcripts",
        font=my_font,
    ).place(relx=0.02, rely=0.86, anchor="w")
    # Define a Checkbox
    case_filter_TE = customtkinter.CTkRadioButton(
        master=sidebar_frame_right, text="False", variable=var_filter_genic, value=1
    ).place(relx=0.02, rely=0.9, anchor="w")
    case_no_filter_TEs = customtkinter.CTkRadioButton(
        master=sidebar_frame_right, text="True", variable=var_filter_genic, value=2
    ).place(relx=0.18, rely=0.9, anchor="w")
    return var_filter_genic


###########################
# window bottom left : Step 3 optional parameters
###########################


def graphical_recip_best_hits_strat1(
    sidebar_frame_south: customtkinter.CTkFrame, my_font: tuple
) -> IntVar:
    """
    Decide whether to use best reciprocal hits or simple reciprocal hits for annotated genes orthology.

    This function creates radio buttons for the user to choose whether to use reciprocal best hits
    (True) or simple reciprocal hits (False) for orthology determination in annotated genes. The user's
    choice is stored in an `IntVar` instance.

    Args:
        sidebar_frame_south (customtkinter.CTkFrame): The GUI frame where the radio buttons are displayed.
        my_font (tuple): The font configuration for the label text.

    Returns:
        IntVar: A variable holding the user's choice for reciprocal best hits or simple reciprocal hits.
    """
    # Define empty variables
    var_rec_best_hit = IntVar(value=1)
    label = customtkinter.CTkLabel(
        sidebar_frame_south, text="Reciprocal BEST hits", font=my_font
    ).place(relx=0.02, rely=0.1, anchor="w")
    # Define a Checkbox
    case_filter_TE = customtkinter.CTkRadioButton(
        master=sidebar_frame_south, text="False", variable=var_rec_best_hit, value=1
    ).place(relx=0.02, rely=0.2, anchor="w")
    case_no_filter_TEs = customtkinter.CTkRadioButton(
        master=sidebar_frame_south, text="True", variable=var_rec_best_hit, value=2
    ).place(relx=0.18, rely=0.2, anchor="w")
    return var_rec_best_hit


def graphical_synteny_window_strat1(
    sidebar_frame_south: customtkinter.CTkFrame, my_font: tuple
) -> IntVar:
    """
    Allows the user to choose the synteny window size.

    This function provides a Spinbox for the user to select a synteny window size (from 0 to 5).
    The value is stored in an `IntVar` which can later be accessed for further processing.

    Args:
        sidebar_frame_south (customtkinter.CTkFrame): The GUI frame where the Spinbox is placed.
        my_font (tuple): The font configuration for the label text.

    Returns:
        IntVar: A variable holding the synteny window size selected by the user.
    """
    label = customtkinter.CTkLabel(
        sidebar_frame_south, text="Synteny window (No synteny : 0)", font=my_font
    ).place(relx=0.02, rely=0.36, anchor="w")
    var_synteny = IntVar(value=2)
    choice_TPM = Spinbox(
        sidebar_frame_south,
        from_=0,
        to=5,
        increment=1,
        textvariable=var_synteny,
        wrap=True,
    ).place(relx=0.02, rely=0.46, anchor="w")
    return var_synteny


def graphical_slider_event_strat1(
    value: float, window: Tk, dico: dict
) -> None:  ###### the value was already in the window
    """
    Handles the event when the slider is moved and updates the dictionary with the selected value.

    This function updates the dictionary `dico` with the value selected from the slider,
    and displays the selected value as a label on the provided window.

    Args:
        value (float): The value from the slider, typically a float, which is converted to an integer.
        window (Tk): The Tkinter window where the label displaying the selected value will be placed.
        dico (Dict[str, int]): The dictionary that holds various configurations, which will be updated with the premature stop value.

    Returns:
        None
    """
    dico["premature_stop"] = int(value)
    label = customtkinter.CTkLabel(
        window,
        text=str(int(value)),
        width=30,
        height=25,
        text_color="#37589c",
        corner_radius=4,
        fg_color="white",
        font=("Arial", 15),
    ).place(relx=0.5, rely=0.68)


def graphical_percent_pos_stop_strat1(
    sidebar_frame_south: customtkinter.CTkFrame, my_font: tuple, dico_file_and_dir: dict
) -> None:
    """
    Allows the user to choose the percentage of a homologous hit where a premature stop codon is validated.

    This function creates a slider to adjust the percentage threshold for identifying premature stop codons.
    It also displays the current percentage as a label that updates when the slider is moved.

    Args:
        sidebar_frame_south (Tk): The Tkinter frame where the slider and labels are placed.
        my_font (str): The font used for text elements in the window.
        dico_file_and_dir (Dict[str, int]): A dictionary where the premature stop codon percentage value will be stored.

    Returns:
        None
    """
    label = customtkinter.CTkLabel(
        sidebar_frame_south, text="Premature STOP codons (%)", font=my_font
    ).place(relx=0.02, rely=0.63, anchor="w")
    label = customtkinter.CTkLabel(
        sidebar_frame_south,
        text="",
        width=30,
        height=25,
        text_color="#37589c",
        corner_radius=4,
        fg_color="white",
        font=("Arial", 15),
    ).place(relx=0.5, rely=0.68)
    label = customtkinter.CTkLabel(
        sidebar_frame_south,
        text="50",
        width=30,
        height=25,
        text_color="#37589c",
        corner_radius=4,
        fg_color="white",
        font=("Arial", 15),
    ).place(relx=0.5, rely=0.68)
    slider = customtkinter.CTkSlider(
        master=sidebar_frame_south,
        from_=0,
        to=100,
        command=lambda value,
        window=sidebar_frame_south,
        dico=dico_file_and_dir: graphical_slider_event_strat1(value, window, dico),
    ).place(relx=0.02, rely=0.73, anchor="w")


###########################
# Main
###########################


def my_graphical_interface_strategy1() -> dict:
    """
    Creates and displays the graphical user interface (GUI) for strategy 1 using the customtkinter library.

    This function initializes a GUI window, sets up frames, buttons, and input fields,
    and allows users to configure various parameters related to transcript analysis.
    It stores user selections in a dictionary and returns it upon closing the interface.

    Returns:
        dict: A dictionary containing user-defined options and settings, including:
            - strategy (str): The strategy identifier ("1").
            - query (str): User-defined query.
            - path_to_genome_repository (str): Path to the genome repository.
            - path_to_transcriptome_repository (str): Path to the transcriptome repository.
            - link_database_outgroup_prot (str): Link to the outgroup protein database.
            - link_database_outgroup_nucl (str): Link to the outgroup nucleotide database.
            - TPM_threeshold (float): Threshold for transcript per million (TPM).
            - transcript_overlap (list): List of transcript overlap types (e.g., "intergenic", "intronic").
            - ORFs_choice (list): List defining ORF selection preferences.
            - filter_genic (bool): Whether to filter genic transcripts.
            - filter_TE (str): Whether to filter transposable elements.
            - rm_undir_transc (str): Whether to remove undirected transcripts.
            - parameters_database_prot (dict): Dictionary containing protein database search parameters.
            - parameters_database_nucl (dict): Dictionary containing nucleotide database search parameters.
            - rec_best_hit (str): Whether to use reciprocal best-hit filtering.
            - synteny_window (int): Size of the synteny window.
            - premature_stop (int): Value for premature stop codon handling.

    Note:
        The function relies on various helper functions (e.g., `graphical_enter_mandatory_strat1`,
        `graphical_TPM_threeshold_strat1`, `graphical_synteny_window_strat1`) to populate the interface
        and retrieve user inputs.

    """
    win = customtkinter.CTk()
    win.title("DESwoMAN")

    img = PhotoImage(file="DESWOMAN_logo_small_transp.png")
    canvas = Canvas(win)
    canvas.configure(
        width=img.width(),
        height=img.height(),
        bg="#222323",
        highlightbackground="#252627",
    )  #
    canvas.create_image(img.width() / 2, img.height() / 2, image=img)
    canvas.place(relx=0.2, rely=0.9, anchor=CENTER)

    customtkinter.set_appearance_mode("dark")
    my_font = customtkinter.CTkFont(
        family="Courier", size=15, weight="bold", underline=False, overstrike=False
    )

    # Set the geometry of Tkinter frame
    win.geometry("1000x690")

    sidebar_frame_left = customtkinter.CTkFrame(
        master=win, width=450, height=280, corner_radius=5
    )
    sidebar_frame_left.place(relx=0.03, rely=0.21, anchor="w")
    sidebar_frame_left.grid_rowconfigure(4, weight=1)

    sidebar_frame_right = customtkinter.CTkFrame(
        master=win, width=420, height=600, corner_radius=5
    )
    sidebar_frame_right.place(relx=0.54, rely=0.445, anchor="w")
    sidebar_frame_right.grid_rowconfigure(4, weight=1)

    sidebar_frame_south = customtkinter.CTkFrame(
        master=win, width=450, height=240, corner_radius=5
    )
    sidebar_frame_south.place(relx=0.03, rely=0.65, anchor="w")
    sidebar_frame_south.grid_rowconfigure(4, weight=1)

    closed_start = customtkinter.CTkFrame(
        master=win, width=100, height=50, corner_radius=5, fg_color="#789bbd"
    )
    closed_start.place(relx=0.05, rely=0.91, anchor="w")
    closed_start.grid_rowconfigure(4, weight=1)
    label = customtkinter.CTkLabel(
        closed_start,
        text="RUN",
        width=30,
        height=25,
        corner_radius=4,
        font=("Courier", 16),
    ).place(relx=0.32, rely=0.3)

    radiobutton_frame = customtkinter.CTkFrame(
        win, width=440, height=80, corner_radius=5, fg_color="#484747"
    )
    radiobutton_frame.place(
        relx=0.035,
        rely=0.08,
        anchor="w",
    )

    ### variables

    dico_file_and_dir = {
        "strategy": "1",
        "query": "",
        "path_to_genome_repository": "",
        "path_to_transcriptome_repository": "",
        "link_database_outgroup_prot": "",
        "link_database_outgroup_nucl": "",
        "TPM_threeshold": 0.5,
        "transcript_overlap": ["intergenic"],
        "ORFs_choice": [["longest"], ["duplicate_handle"], ["utr_size", 0, 0]],
        "filter_genic": False,
        "filter_TE": "False",
        "rm_undir_transc": "False",
        "parameters_database_prot": {"type": "blastp", "mode": "--more-sensitive"},
        "parameters_database_nucl": {
            "type": "blastn",
            "e_value": "0.01",
            "coverage": None,
            "strand": None,
        },
        "rec_best_hit": "False",
        "synteny_window": 2,
        "premature_stop": 50,
    }

    ## function call

    # databases
    value_ref_name = graphical_enter_mandatory_strat1(
        win, radiobutton_frame, my_font, dico_file_and_dir
    )
    var_blastp_mode, var_blastnucl_mode = graphical_enter_databases_strat1(
        sidebar_frame_left, my_font, dico_file_and_dir
    )

    # transcripts options
    var_TPM = graphical_TPM_threeshold_strat1(sidebar_frame_right, my_font)
    var_filter_TE = graphical_filter_TE_strat1(sidebar_frame_right, my_font)
    var_undir_transcr = graphical_remove_undir_transcripts_strat1(
        sidebar_frame_right, my_font
    )
    var_intergenic, var_intronic, var_antisense, var_genic = (
        graphical_transcript_location_strat1(sidebar_frame_right, my_font)
    )
    var_orf_choice = graphical_ORFs_choice_strat1(sidebar_frame_right, my_font)
    var_5UTR = graphical_min_size_5UTR_strat1(sidebar_frame_right, my_font)
    var_3UTR = graphical_min_size_3UTR_strat1(sidebar_frame_right, my_font)
    var_filter_genic = graphical_filter_genic_transcripts_strat1(
        sidebar_frame_right, my_font
    )
    graphical_percent_pos_stop_strat1(sidebar_frame_south, my_font, dico_file_and_dir)

    # step 3
    var_rec_best_hit = graphical_recip_best_hits_strat1(sidebar_frame_south, my_font)
    var_synteny = graphical_synteny_window_strat1(sidebar_frame_south, my_font)

    ### outputs
    win.mainloop()
    dico_file_and_dir["query"] = value_ref_name.get()
    dico_file_and_dir["TPM_threeshold"] = float(var_TPM.get())
    if (
        var_intergenic.get() == 0
        and "intergenic" in dico_file_and_dir["transcript_overlap"]
    ):
        dico_file_and_dir["transcript_overlap"].remove("intergenic")

    if (
        var_intergenic.get() == 1
        and "intergenic" not in dico_file_and_dir["transcript_overlap"]
    ):
        dico_file_and_dir["transcript_overlap"].append("intergenic")

    if (
        var_intronic.get() == 0
        and "intronic" in dico_file_and_dir["transcript_overlap"]
    ):
        dico_file_and_dir["transcript_overlap"].remove("intronic")

    if (
        var_intronic.get() == 1
        and "intronic" not in dico_file_and_dir["transcript_overlap"]
    ):
        dico_file_and_dir["transcript_overlap"].append("intronic")

    if (
        var_antisense.get() == 0
        and "antisense" in dico_file_and_dir["transcript_overlap"]
    ):
        dico_file_and_dir["transcript_overlap"].remove("antisense")

    if (
        var_antisense.get() == 1
        and "antisense" not in dico_file_and_dir["transcript_overlap"]
    ):
        dico_file_and_dir["transcript_overlap"].append("antisense")

    if var_genic.get() == 0 and "genic" in dico_file_and_dir["transcript_overlap"]:
        dico_file_and_dir["transcript_overlap"].remove("genic")

    if var_genic.get() == 1 and "genic" not in dico_file_and_dir["transcript_overlap"]:
        dico_file_and_dir["transcript_overlap"].append("genic")
    dico_file_and_dir["TPM_threeshold"] = float(var_TPM.get())

    if var_orf_choice.get() == 1 and ["all"] not in dico_file_and_dir["ORFs_choice"]:
        dico_file_and_dir["ORFs_choice"].append(["all"])
        if ["longest"] in dico_file_and_dir["ORFs_choice"]:
            dico_file_and_dir["ORFs_choice"].remove(["longest"])
        if ["kozac_highest"] in dico_file_and_dir["ORFs_choice"]:
            dico_file_and_dir["ORFs_choice"].remove(["kozac_highest"])
        if ["start_first"] in dico_file_and_dir["ORFs_choice"]:
            dico_file_and_dir["ORFs_choice"].remove(["start_first"])

    if (
        var_orf_choice.get() == 2
        and ["longest"] not in dico_file_and_dir["ORFs_choice"]
    ):
        dico_file_and_dir["ORFs_choice"].append(["longest"])
        if ["all"] in dico_file_and_dir["ORFs_choice"]:
            dico_file_and_dir["ORFs_choice"].remove(["all"])
        if ["kozac_highest"] in dico_file_and_dir["ORFs_choice"]:
            dico_file_and_dir["ORFs_choice"].remove(["kozac_highest"])
        if ["start_first"] in dico_file_and_dir["ORFs_choice"]:
            dico_file_and_dir["ORFs_choice"].remove(["start_first"])

    if (
        var_orf_choice.get() == 3
        and ["kozac_highest"] not in dico_file_and_dir["ORFs_choice"]
    ):
        dico_file_and_dir["ORFs_choice"].append(["kozac_highest"])
        if ["all"] in dico_file_and_dir["ORFs_choice"]:
            dico_file_and_dir["ORFs_choice"].remove(["all"])
        if ["longest"] in dico_file_and_dir["ORFs_choice"]:
            dico_file_and_dir["ORFs_choice"].remove(["longest"])
        if ["start_first"] in dico_file_and_dir["ORFs_choice"]:
            dico_file_and_dir["ORFs_choice"].remove(["start_first"])

    if (
        var_orf_choice.get() == 4
        and ["start_first"] not in dico_file_and_dir["ORFs_choice"]
    ):
        dico_file_and_dir["ORFs_choice"].append(["start_first"])
        if ["all"] in dico_file_and_dir["ORFs_choice"]:
            dico_file_and_dir["ORFs_choice"].remove(["all"])
        if ["longest"] in dico_file_and_dir["ORFs_choice"]:
            dico_file_and_dir["ORFs_choice"].remove(["longest"])
        if ["kozac_highest"] in dico_file_and_dir["ORFs_choice"]:
            dico_file_and_dir["ORFs_choice"].remove(["kozac_highest"])

    if var_filter_genic.get() == 1:
        dico_file_and_dir["filter_genic"] = False
    else:
        dico_file_and_dir["filter_genic"] = True

    if var_filter_TE.get() == 1:
        dico_file_and_dir["filter_TE"] = "False"
    else:
        dico_file_and_dir["filter_TE"] = "True"

    if var_undir_transcr.get() == 1:
        dico_file_and_dir["rm_undir_transc"] = "False"
    else:
        dico_file_and_dir["rm_undir_transc"] = "True"

    if var_5UTR.get() != 0:
        already_present = False
        for i in range(len(dico_file_and_dir["ORFs_choice"])):
            if "utr_size" in dico_file_and_dir["ORFs_choice"][i]:
                dico_file_and_dir["ORFs_choice"][i][1] = var_5UTR.get()
                already_present = True
        if already_present == False:
            utr_5_value = var_5UTR.get()
            dico_file_and_dir["ORFs_choice"].append(["utr_size", utr_5_value, 0])

    if var_3UTR.get() != 0:
        already_present = False
        for i in range(len(dico_file_and_dir["ORFs_choice"])):
            if "utr_size" in dico_file_and_dir["ORFs_choice"][i]:
                dico_file_and_dir["ORFs_choice"][i][2] = var_3UTR.get()
                already_present = True
        if already_present == False:
            utr_3_value = var_3UTR.get()
            dico_file_and_dir["ORFs_choice"].append(["utr_size", 0, utr_3_value])

    if var_blastp_mode.get() == 1:
        dico_file_and_dir["parameters_database_prot"]["mode"] = "--sensitive"
    elif var_blastp_mode.get() == 2:
        dico_file_and_dir["parameters_database_prot"]["mode"] = "--more-sensitive"
    elif var_blastp_mode.get() == 3:
        dico_file_and_dir["parameters_database_prot"]["mode"] = "--very-sensitive"
    else:
        dico_file_and_dir["parameters_database_prot"]["mode"] = "--ultra-sensitive"

    if var_blastnucl_mode.get() == 1:
        dico_file_and_dir["parameters_database_nucl"]["e_value"] = "0.1"
    elif var_blastnucl_mode.get() == 2:
        dico_file_and_dir["parameters_database_nucl"]["e_value"] = "0.01"
    elif var_blastnucl_mode.get() == 3:
        dico_file_and_dir["parameters_database_nucl"]["e_value"] = "0.001"
    else:
        dico_file_and_dir["parameters_database_nucl"]["e_value"] = "0.0001"

    if var_rec_best_hit.get() == 1:
        dico_file_and_dir["rec_best_hit"] = "False"
    else:
        dico_file_and_dir["rec_best_hit"] = "True"

    dico_file_and_dir["synteny_window"] = var_synteny.get()
    return dico_file_and_dir
