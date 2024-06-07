import shared_data as sd


pdb_lnc = f"{sd.SYSTEMS['LNC']}/antibiotic_aligned.pdb"
for system_name, system_path in sd.SYSTEMS.items():
    pdb_in = f"{system_path}/common_atoms.pdb"
    pdb_out = f"{system_path}/added_antibiotic.pdb"
    with open(pdb_in, "r") as file_in, open(pdb_lnc, "r") as file_lnc, open(pdb_out, "w") as file_out:
        lines_in = file_in.readlines()
        lines_lnc = file_lnc.readlines()
        index = int(lines_in[-3][6:11]) + 1
        for line in lines_in:
            if "END" not in line and "TER" not in line:
                file_out.write(line)
        for line in lines_lnc:
            if "ATOM" in line and line[21].strip() == "9":
                line1 = line[:6]
                line2 = line[11:]
                index_string = str(index).rjust(5)
                line = line1 + index_string + line2
                file_out.write(line)
                index += 1
        file_out.write("END")