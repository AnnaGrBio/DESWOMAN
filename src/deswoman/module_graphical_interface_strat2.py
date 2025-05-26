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
# window left : Input data
###########################


def graphical_select_genome_folder_strat2(
    dico_file_and_dir: dict, radiobutton_frame: Tk, win: Tk, value_ref_name: StringVar
) -> None:
    """
    Opens a file dialog to select a genome folder and updates the provided dictionary.

    This function allows the user to select a genome repository folder using a file dialog.
    Once selected, it updates `dico_file_and_dir` with the folder path and visually confirms the selection.
    If both `path_to_transcriptome_repository` and a reference name are provided, it enables the "RUN" button.

    Args:
        dico_file_and_dir (Dict[str, str]): Dictionary storing file paths and directories.
        radiobutton_frame (Tk): The Tkinter frame where the confirmation label is placed.
        win (Tk): The main window where the "RUN" button is placed.
        value_ref_name (customtkinter.StringVar): Variable holding the reference name.

    Returns:
        None
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
    ).place(relx=0.38, rely=0.65)
    name_ref = value_ref_name.get()
    if dico_file_and_dir["path_to_transcriptome_repository"] != "" and name_ref != "":
        button_DESMAN_run = customtkinter.CTkButton(
            master=win,
            width=100,
            height=50,
            text="RUN",
            command=lambda: graphical_Close_strat2(win),
        ).place(relx=0.05, rely=0.81, anchor="w")


def graphical_select_transcriptome_folder_strat2(
    dico_file_and_dir: dict, radiobutton_frame: Tk, win: Tk, value_ref_name: StringVar
) -> None:
    """
    Opens a file dialog to select a transcriptome folder and updates the provided dictionary.

    Once a folder is selected, the function updates `dico_file_and_dir` with the folder path
    and visually confirms the selection. If both `path_to_genome_repository` and a reference
    name are provided, the "RUN" button is enabled.

    Args:
        dico_file_and_dir (Dict[str, str]): Dictionary storing file paths and directories.
        radiobutton_frame (Tk): The Tkinter frame where the confirmation label is placed.
        win (Tk): The main window where the "RUN" button is placed.
        value_ref_name (customtkinter.StringVar): Variable holding the reference name.

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
    ).place(relx=0.905, rely=0.65)
    name_ref = value_ref_name.get()
    if dico_file_and_dir["path_to_genome_repository"] != "" and name_ref != "":
        button_DESMAN_run = customtkinter.CTkButton(
            master=win,
            width=100,
            height=50,
            text="RUN",
            command=lambda: graphical_Close_strat2(win),
        ).place(relx=0.05, rely=0.81, anchor="w")


def graphical_Close_strat2(win: Tk) -> None:
    """
    Closes the given Tkinter window.

    Args:
        win (Tk): The Tkinter window to be closed.

    Returns:
        None
    """
    win.destroy()


def graphical_activate_with_text_strat2(win: Tk, dico_file_and_dir: dict) -> None:
    """
    Activates the RUN button if genome and transcriptome paths are selected.

    Args:
        win: The Tkinter window where the button should be displayed.
        dico_file_and_dir (dict): Dictionary containing paths to genome and transcriptome repositories.

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
            command=lambda: graphical_Close_strat2(win),
        ).place(relx=0.05, rely=0.81, anchor="w")


def graphical_enter_mandatory_strat2(
    win: Tk, radiobutton_frame: Tk, my_font: tuple, dico_file_and_dir: dict
) -> StringVar:
    """
    Creates input fields and buttons for selecting genome and transcriptome folders,
    as well as entering a query name. Also updates the UI dynamically based on user input.

    Args:
        win (Tk): The main Tkinter window.
        radiobutton_frame (Frame): The frame containing the UI elements.
        my_font (tuple): The font settings for UI elements.
        dico_file_and_dir (dict): A dictionary to store selected file paths.

    Returns:
        StringVar: A Tkinter variable containing the query name input by the user.
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
        lambda e: graphical_activate_with_text_strat2(win, dico_file_and_dir),
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
    ).place(relx=0.38, rely=0.65)
    label = customtkinter.CTkLabel(
        radiobutton_frame,
        text="",
        font=my_font,
        corner_radius=4,
        width=30,
        height=25,
        fg_color="white",
    ).place(relx=0.905, rely=0.65)
    button_genome_dir = customtkinter.CTkButton(
        master=radiobutton_frame,
        width=150,
        height=25,
        text="Genome directory",
        command=lambda: graphical_select_genome_folder_strat2(
            dico_file_and_dir, radiobutton_frame, win, value_ref_name
        ),
    ).place(relx=0.02, rely=0.75, anchor="w")
    button_transcriptome_dir = customtkinter.CTkButton(
        master=radiobutton_frame,
        width=150,
        height=25,
        text="Transcriptome directory",
        command=lambda: graphical_select_transcriptome_folder_strat2(
            dico_file_and_dir, radiobutton_frame, win, value_ref_name
        ),
    ).place(relx=0.54, rely=0.75, anchor="w")
    return value_ref_name


###########################
### dataset for Step 2 Blast homology search
###########################


def graphical_select_prot_database_strat2(
    radiobutton_frame2: Tk, dico_file_and_dir: dict
) -> None:
    """
    Opens a file dialog to select a protein dataset (FASTA format) for homology search in Step 2.
    Updates the provided dictionary with the selected file path and displays a confirmation label.

    Args:
        radiobutton_frame2 (Frame): The frame where the UI elements are displayed.
        dico_file_and_dir (dict): A dictionary to store the selected protein database path.

    Returns:
        None
    """
    filepath = askopenfilename(
        title="Select a dataset",
        filetypes=[("fasta file", ".fa .fna .fasta"), ("all files", ".*")],
    )
    dico_file_and_dir["link_database_outgroup_prot"] = filepath
    label = customtkinter.CTkLabel(
        radiobutton_frame2,
        text="V",
        width=30,
        height=25,
        text_color="#37589c",
        corner_radius=4,
        font=("Arial", 15),
    ).place(relx=0.36, rely=0.418)


def graphical_select_nucl_database_strat2(
    sidebar_frame_left: Tk, dico_file_and_dir: dict
) -> None:
    """
    Opens a file dialog to select a nucleotide dataset (FASTA format) for homology search in Step 2.
    Updates the provided dictionary with the selected file path and displays a confirmation label.

    Args:
        sidebar_frame_left (Frame): The frame where the UI elements are displayed.
        dico_file_and_dir (dict): A dictionary to store the selected nucleotide database path.

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


def graphical_enter_databases_strat2(
    sidebar_frame_left: Tk, my_font: tuple, dico_file_and_dir: dict
) -> tuple:
    """
    Creates a graphical interface for selecting BLAST protein and nucleotide databases
    and configuring sensitivity levels for both BLASTp and BLASTn searches.

    Args:
        sidebar_frame_left (Frame): The frame where the UI elements are displayed.
        my_font (tuple): The font style used for labels and buttons.
        dico_file_and_dir (dict): A dictionary to store selected database paths.

    Returns:
        tuple[IntVar, IntVar]: Variables representing sensitivity levels for BLASTp and BLASTn searches.
    """
    button_prot_file = customtkinter.CTkButton(
        master=sidebar_frame_left,
        width=150,
        height=25,
        text="BLAST protein database",
        command=lambda: graphical_select_prot_database_strat2(
            sidebar_frame_left, dico_file_and_dir
        ),
    ).place(relx=0.02, rely=0.18, anchor="w")
    button_nucl_file = customtkinter.CTkButton(
        master=sidebar_frame_left,
        width=150,
        height=25,
        text="BLAST nucl database",
        command=lambda: graphical_select_nucl_database_strat2(
            sidebar_frame_left, dico_file_and_dir
        ),
    ).place(relx=0.02, rely=0.68, anchor="w")
    var_blastp_mode = IntVar(value=2)
    # Define a Checkbox
    case_1 = customtkinter.CTkRadioButton(
        master=sidebar_frame_left, text="sensitive", variable=var_blastp_mode, value=1
    ).place(relx=0.02, rely=0.38, anchor="w")
    case_2 = customtkinter.CTkRadioButton(
        master=sidebar_frame_left,
        text="more sensitive",
        variable=var_blastp_mode,
        value=2,
    ).place(relx=0.22, rely=0.38, anchor="w")
    case_3 = customtkinter.CTkRadioButton(
        master=sidebar_frame_left,
        text="very sensitive",
        variable=var_blastp_mode,
        value=3,
    ).place(relx=0.49, rely=0.38, anchor="w")
    case_4 = customtkinter.CTkRadioButton(
        master=sidebar_frame_left,
        text="ultra sensitive",
        variable=var_blastp_mode,
        value=4,
    ).place(relx=0.75, rely=0.38, anchor="w")

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


def graphical_TPM_threeshold_strat2(
    sidebar_frame_right: Tk, my_font: tuple
) -> StringVar:
    """
    Creates a graphical interface for defining the TPM (Transcripts Per Million) threshold.

    Args:
        sidebar_frame_right (Frame): The frame where the UI elements are displayed.
        my_font (tuple): The font style used for labels.

    Returns:
        StringVar: A variable storing the selected TPM threshold.
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


def graphical_filter_TE_strat2(sidebar_frame_right: Tk, my_font: tuple) -> IntVar:
    """
    Creates a graphical interface for deciding whether to filter Transposable Elements (TEs).

    Args:
        sidebar_frame_right (Frame): The frame where the UI elements will be placed.
        my_font (tuple): The font style used for labels.

    Returns:
        customtkinter.IntVar: A variable storing the choice of whether to filter TEs (1 for False, 2 for True).
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


def graphical_remove_undir_transcripts_strat2(
    sidebar_frame_right: Tk, my_font: tuple
) -> IntVar:
    """
    Creates a graphical interface for deciding whether to remove transcripts with an undefined (".") direction.

    Args:
        sidebar_frame_right (Frame): The frame where the UI elements will be placed.
        my_font (tuple): The font style used for the labels.

    Returns:
        customtkinter.IntVar: A variable that stores the choice to remove or keep unoriented transcripts (1 for False, 2 for True).
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


def graphical_transcript_location_strat2(sidebar_frame_right, my_font):
    """
    Decide which transcripts to keep according to their genomic location
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


def graphical_ORFs_choice_strat2(sidebar_frame_right: Tk, my_font: tuple) -> tuple:
    """
    Creates a graphical interface to decide which transcripts to keep based on their genomic location.

    Args:
        sidebar_frame_right (Frame): The frame where the UI elements will be placed.
        my_font (tuple): The font style used for the labels.

    Returns:
        tuple: A tuple of IntVar variables representing whether the user wants to keep the transcripts in the specified locations.
               (intergenic, intronic, antisense, genic)
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


def graphical_min_size_5UTR_strat2(
    sidebar_frame_right: Tk, my_font: tuple
) -> StringVar:
    """
    Creates a graphical interface to define the minimum length of the 5'UTR.

    Args:
        sidebar_frame_right (customtkinter.CTkFrame): The frame where the UI elements will be placed.
        my_font (tuple): The font style used for the label.

    Returns:
        StringVar: A StringVar object containing the minimum 5'UTR size selected by the user.
    """
    label = customtkinter.CTkLabel(
        sidebar_frame_right, text="Min 5'UTR size (nucl)", font=my_font
    ).place(relx=0.02, rely=0.75, anchor="w")
    var_5UTR = StringVar(value=0)
    choice_5prime_min_size = Spinbox(
        sidebar_frame_right, from_=0, to=100, textvariable=var_5UTR, wrap=True
    ).place(relx=0.02, rely=0.79, anchor="w")
    return var_5UTR


def graphical_min_size_3UTR_strat2(
    sidebar_frame_right: Tk, my_font: tuple
) -> StringVar:
    """
    Creates a graphical interface to define the minimum length of the 3'UTR.

    Args:
        sidebar_frame_right (customtkinter.CTkFrame): The frame where the UI elements will be placed.
        my_font (tuple): The font style used for the label.

    Returns:
        StringVar: A StringVar object containing the minimum 3'UTR size selected by the user.
    """
    label = customtkinter.CTkLabel(
        sidebar_frame_right, text="Min 3'UTR size (nucl)", font=my_font
    ).place(relx=0.54, rely=0.75, anchor="w")
    var_3UTR = StringVar(value=0)
    choice_3prime_min_size = Spinbox(
        sidebar_frame_right, from_=0, to=100, textvariable=var_3UTR, wrap=True
    ).place(relx=0.56, rely=0.79, anchor="w")
    return var_3UTR


def graphical_filter_genic_transcripts_strat2(
    sidebar_frame_right: Tk, my_font: tuple
) -> IntVar:
    """
    Decide whether to filter transcripts that correspond to new splicing events of known genes,
    even though they do not overlap with the gene.

    Args:
        sidebar_frame_right (customtkinter.CTkFrame): The frame where the UI elements will be placed.
        my_font (tuple): The font style used for the label.

    Returns:
        IntVar: An IntVar object representing whether to filter the splice variants of genic transcripts.
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


def my_graphical_interface_strategy2():
    """
    Creates and displays the graphical user interface (GUI) for strategy 2 using the customtkinter library.

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

    Note:
        The function relies on various helper functions (e.g., `graphical_enter_mandatory_strat1`,
        `graphical_TPM_threeshold_strat1`, `graphical_synteny_window_strat1`) to populate the interface
        and retrieve user inputs.

    """
    # Create an instance of tkinter frame
    win = customtkinter.CTk()
    win.title("DESwoMAN")

    img = PhotoImage(file="DESWOMAN_logo_small_transp.png")  #
    canvas = Canvas(win)
    canvas.configure(
        width=img.width(),
        height=img.height(),
        bg="#222323",
        highlightbackground="#252627",
    )
    canvas.create_image(img.width() / 2, img.height() / 2, image=img)
    canvas.place(relx=0.2, rely=0.8, anchor=CENTER)

    customtkinter.set_appearance_mode("dark")
    my_font = customtkinter.CTkFont(
        family="Courier", size=15, weight="bold", underline=False, overstrike=False
    )

    win.geometry("1000x590")

    sidebar_frame_left = customtkinter.CTkFrame(
        master=win, width=480, height=360, corner_radius=5
    )
    sidebar_frame_left.place(relx=0.03, rely=0.32, anchor="w")
    sidebar_frame_left.grid_rowconfigure(4, weight=1)

    sidebar_frame_right = customtkinter.CTkFrame(
        master=win, width=420, height=570, corner_radius=5
    )
    sidebar_frame_right.place(relx=0.54, rely=0.5, anchor="w")
    sidebar_frame_right.grid_rowconfigure(4, weight=1)

    closed_start = customtkinter.CTkFrame(
        master=win, width=100, height=50, corner_radius=5, fg_color="#789bbd"
    )
    closed_start.place(relx=0.05, rely=0.81, anchor="w")
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
        win, width=450, height=120, corner_radius=5, fg_color="#484747"
    )
    radiobutton_frame.place(
        relx=0.045,
        rely=0.14,
        anchor="w",
    )

    radiobutton_frame2 = customtkinter.CTkFrame(
        win, width=450, height=200, corner_radius=5, fg_color="#484747"
    )
    radiobutton_frame2.place(
        relx=0.045,
        rely=0.44,
        anchor="w",
    )

    ### variables

    dico_file_and_dir = {
        "strategy": "2",
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
    }

    ## function call

    # databases
    value_ref_name = graphical_enter_mandatory_strat2(
        win, radiobutton_frame, my_font, dico_file_and_dir
    )
    var_blastp_mode, var_blastnucl_mode = graphical_enter_databases_strat2(
        radiobutton_frame2, my_font, dico_file_and_dir
    )

    # transcripts options
    var_TPM = graphical_TPM_threeshold_strat2(sidebar_frame_right, my_font)
    var_filter_TE = graphical_filter_TE_strat2(sidebar_frame_right, my_font)
    var_undir_transcr = graphical_remove_undir_transcripts_strat2(
        sidebar_frame_right, my_font
    )
    var_intergenic, var_intronic, var_antisense, var_genic = (
        graphical_transcript_location_strat2(sidebar_frame_right, my_font)
    )
    var_orf_choice = graphical_ORFs_choice_strat2(sidebar_frame_right, my_font)
    var_5UTR = graphical_min_size_5UTR_strat2(sidebar_frame_right, my_font)
    var_3UTR = graphical_min_size_3UTR_strat2(sidebar_frame_right, my_font)
    var_filter_genic = graphical_filter_genic_transcripts_strat2(
        sidebar_frame_right, my_font
    )

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
    # print (remove_undir_tr)

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

    return dico_file_and_dir
