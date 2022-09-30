def generate_gene_type():
    name_type = {}

    base_file = open("my_gene_type.txt")
    for line in base_file.readlines():
        # info = line.split()[6:]
        info = line.split()
        g_name = info[0]
        g_type = info[1]
        g_name = g_name.lower()
        if g_name not in name_type:
            if "pseudogene" in g_type:
                g_type = "pseudogene"
            name_type[g_name] = g_type
    base_file.close()

    # base_file_2 = open("combined.gene.bed+2")
    # for line in base_file_2.readlines():
    #     # info = line.split()[6:]
    #     info = line.split()
    #     g_name = info[-2]
    #     g_type = info[-1]
    #     g_name = g_name.lower()
    #     if g_name not in name_type:
    #         if "pseudogene" in g_type:
    #             g_type = "pseudogene"
    #         name_type[g_name] = g_type
    # base_file_2.close()

    gene_type_file = open("searching_results.txt")
    # sample content:
    #   Searching: loc102723975 ...		type:  unknown_type
    for line in gene_type_file.readlines():
        info = line.split()
        g_name = info[1]
        g_type = info[4]
        if g_name not in name_type:
            name_type[g_name] = g_type
    gene_type_file.close()

    my_gene_type = open("my_gene_type_2.txt", "w")
    for name in name_type.keys():
        g_type = name_type[name]
        line = "%s %s\n" % (name, g_type)
        my_gene_type.write(line)
    my_gene_type.close()

    print("Successfully generate to => my_gene_type.txt !")


if __name__ == '__main__':
    generate_gene_type()
