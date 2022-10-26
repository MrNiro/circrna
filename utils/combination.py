import os

from selenium import webdriver
from selenium.common.exceptions import TimeoutException
from selenium.webdriver.common.by import By
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.support.wait import WebDriverWait

chrome_option = webdriver.ChromeOptions()
chrome_option.add_argument(argument='--headless')
browser = webdriver.Chrome(options=chrome_option)   # options=chrome_option
wait = WebDriverWait(browser, 2)  # 最长等待时间为3S
# wait_v = WebDriverWait(browser, 20)
browser.maximize_window()


def crawler(g_name):
    if g_name == "na":
        return "unknown_type"

    base_url = "https://www.genecards.org/cgi-bin/carddisp.pl?gene="
    url = base_url + g_name
    browser.get(url)
    print("Searching:", g_name, end=" ...")
    try:
        g_type = wait.until(EC.presence_of_element_located(
            (By.CSS_SELECTOR, "#aliases_descriptions > div.row "
                              "> div.col-xs-12.col-md-9 > div:nth-child(2) > div > div"))).text
    except TimeoutException:
        try:
            temp = wait.until(EC.presence_of_element_located(
                (By.CSS_SELECTOR, "#cardPage > div:nth-child(1) > div > div > div.gc-section-header.gc-first "
                                  "> div.col-xs-8.col-md-10.col-sm-9.gene-symbol-description.row-no-padding "
                                  "> div > span.gc-category"))).text
            if temp == "Protein Coding":
                g_type = "protein_coding"
            elif temp == "Pseudogene":
                g_type = "pseudogene"
            else:
                print("----------------------------------", temp)
                g_type = temp

        except TimeoutException:
            g_type = "unknown_type"
    print("\t\ttype: ", g_type)
    return g_type


def combine_circexplore(circexplore_path):
    all_detected_circrna = {}

    for _, dirs, _ in os.walk(circexplore_path):
        for dir_name in dirs:
            temp_path = circexplore_path + dir_name
            for _, _, filenames in os.walk(temp_path):
                for f in filenames:
                    file_path = temp_path + "/" + f
                    bed_file = open(file_path)
                    for line in bed_file.readlines():
                        info = line.strip().split()
                        # rna_id = info[3]
                        rna_id = "%s:%s-%s:%s" % (info[0], info[1], info[2], info[5])

                        rna_type = info[13]
                        rna_name = info[14]
                        if "," in rna_name:
                            rna_name = rna_name.split(",")[0]

                        if rna_id not in all_detected_circrna:
                            all_detected_circrna[rna_id] = [rna_type, rna_name]
    return all_detected_circrna


def combine_DEA_circrna(all_detected_circrna):
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
                            log2_fold_change = info[2]
                            p_value = info[5]
                            if rna_id not in all_rna:
                                all_rna[rna_id] = [log2_fold_change, p_value]
                            else:
                                all_rna[rna_id].extend([log2_fold_change, p_value])
                        break

    for each in all_rna.keys():
        if each in all_detected_circrna:
            all_rna[each].append(all_detected_circrna[each])
        else:
            print(each)
            print(all_rna[each])
            raise Exception()
        # print(all_rna[each])
    return all_rna, diff_types


def load_gene_type():
    gene_type_file = open("gene_type/my_gene_type.txt")
    gene_type = {}
    for line in gene_type_file.readlines():
        g_name, g_type = line.strip().split(" ")
        gene_type[g_name] = g_type
    gene_type_file.close()

    return gene_type


def load_normalize_count():
    f = open("../results/differential_expression/circRNA/DESeq2_normalized_counts.txt")
    rna_count = {}
    for line in f.readlines()[1:]:
        info = line.strip().split("\t")
        g_id = info[0]
        g_count = list(map(lambda x: round(float(x)), info[1:]))
        g_count.insert(0, sum(g_count))
        g_count = str(g_count)[1:-1]
        rna_count[g_id] = g_count

    return rna_count


def generate_csv(all_rna, diff_types, rna_count):
    dea_groups_num = len(diff_types)

    combine_result = open("../results/combined_differential_expression.csv", "w")

    # generate csv Title
    csv_title = "circrna_ID,circType,geneName,geneType,"
    for each in diff_types:
        csv_title += ("log2FC_%s,P_value_%s," % (each, each))

    csv_title += "Normalized_Count_Sum,"
    for i in range(110, 116):
        csv_title += ("RW" + str(i) + "_Count,")

    csv_title += "\n"
    combine_result.write(csv_title)

    all_gene_type = load_gene_type()
    for each in all_rna.keys():
        value = all_rna[each]

        # g_type = "unknown_type"
        circ_type = value[-1][0]
        g_name = value[-1][1].lower()
        if g_name in all_gene_type:
            g_type = all_gene_type[g_name]
        else:
            g_type = crawler(g_name)

        line = ("%s," * 4) \
               % (each, circ_type, g_name.upper(), g_type)

        # add the log2FC and p value
        for i in range(dea_groups_num * 2):
            line += "%s," % value[i]
        # line += ("%s," * dea_groups_num * 2) % value

        if each in rna_count:
            line += rna_count[each]
        line += "\n"
        combine_result.write(line)

    combine_result.close()


if __name__ == '__main__':
    detected_circrna = combine_circexplore("../results/circexploer-1026/")
    my_rna_count = load_normalize_count()
    my_rna, my_diff = combine_DEA_circrna(detected_circrna)
    print("Circle RNA number: %d" % len(my_rna))

    generate_csv(my_rna, my_diff, my_rna_count)
