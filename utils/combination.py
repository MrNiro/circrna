import csv
import os


def combine_bed():
    all_bed = {}

    bed_path = "../results/circexplore2_output/"
    for _, dirs, _ in os.walk(bed_path):
        for dir_name in dirs:
            temp_path = bed_path + dir_name
            for _, _, filenames in os.walk(temp_path):
                for f in filenames:
                    file_path = temp_path + "/" + f
                    bed_file = open(file_path)
                    for line in bed_file.readlines():
                        info = line.strip().split()
                        rna_id = info[3]
                        rna_type = info[12]
                        rna_name = info[13]
                        if rna_id not in all_bed:
                            all_bed[rna_id] = (rna_type, rna_name)
    return all_bed


def load_circrna(all_bed):
    all_rna = {}
    diff_types = []
    rna_path = "../results/differential_expression/circRNA/"
    for _, dirs, _ in os.walk(rna_path):
        for dir_name in dirs:
            temp_path = rna_path + dir_name
            diff_types.append(dir_name)
            for _, _, filenames in os.walk(temp_path):
                for f in filenames:
                    if "whole_differential_expression.txt" in f:
                        file_path = temp_path + "/" + f
                        rna_file = open(file_path)
                        for line in rna_file.readlines()[1:]:
                            info = line.strip().split()
                            rna_id = info[0]
                            log2_fold_change = float(info[2])
                            p_value = float(info[5])
                            if rna_id not in all_rna:
                                all_rna[rna_id] = [log2_fold_change, p_value]
                            else:
                                all_rna[rna_id].extend([log2_fold_change, p_value])
                        break

    for each in all_rna.keys():
        if each in all_bed:
            all_rna[each].append(all_bed[each])
        print(all_rna[each])
    return all_rna, diff_types


def generate_csv(all_rna, diff_types):
    combine_result = open("../results/combined_differential_expression.csv", "w")
    csv_title = "circrna_ID,"
    for each in diff_types:
        csv_title += ("log2FoldChange_%s,p_value_%s," % (each, each))
    csv_title += "circType,geneName,geneName2,geneName3\n"
    combine_result.write(csv_title)

    for each in all_rna.keys():
        value = all_rna[each]
        line = "%s,%lf,%lf,%f,%f,%f,%f,%s,%s\n" \
               % (each, value[0], value[1], value[2], value[3], value[4], value[5], value[6][0], value[6][1])
        combine_result.write(line)

    combine_result.close()


if __name__ == '__main__':
    my_bed = combine_bed()
    my_rna, my_diff = load_circrna(my_bed)
    print("Circle RNA number: %d" % len(my_rna))

    generate_csv(my_rna, my_diff)
