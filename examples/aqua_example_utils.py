def find_aquagpusph_exe():
    """Find the AQUAgpusph executable
    """
    import platform
    import os
    if platform.system() == "Windows":
        exe_name = "AQUAgpusph.exe"
    else:
        exe_name = "AQUAgpusph"
    script_folder = os.path.dirname(os.path.realpath(__file__))
    candidates = "../aquagpusph", "../../../bin"
    for candidate in candidates:
        exe_path = os.path.join(script_folder, candidate, exe_name)
        if os.path.isfile(exe_path):
            return os.path.realpath(exe_path)
    # Not found, we can just rely on the PATH environment variable
    return exe_name


def configure(replacements, in_path, out_path="./"):
    """Configure the files on the specified folder

    Parameters:
    -----------

    replacements (dict): Dictionary whith the replacements to be done.
    in_path (str): The folder where the templates are located
    out_path (str): The folrder where the files shall be printed
    """
    import os
    # Add the executable path if it was not already there
    if "AQUAGPUSPH_EXE" not in replacements.keys():
        replacements["AQUAGPUSPH_EXE"] = find_aquagpusph_exe()
    for (in_dir, _, fnames) in os.walk(in_path):
        out_dir = os.path.join(out_path, in_dir.replace(in_path, ""))
        for fname in fnames:
            with open(os.path.join(in_dir, fname), "r") as fin, \
                 open(os.path.join(out_dir, fname), "w") as fout:
                txt = fin.read()
                for k, v in replacements.items():
                    txt = txt.replace('{{' + k + '}}', v)
                fout.write(txt)
