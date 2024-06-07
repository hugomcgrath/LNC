from glob import glob
import utility_functions as uf

for logfile_name in sorted(glob("/home/hmcgrat/LNC/LOGS/*")):
    cf = uf.get_config_file_from_file_path(logfile_name, type_of_file="log")
    all_ok = 0
    with open(logfile_name, "r") as logfile:
        lines = logfile.readlines()
        for line in lines:
            if "Time" in line:
                line = line.strip()
                print(f"config ID {cf.config_index}:\t{line}")
                all_ok = 1
    if not all_ok:
        print(f"config ID {cf.config_index}:\tMissing logfile")
        