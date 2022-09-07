import sys
import os


if __name__ == '__main__':
    data_folder = "/home/ubuntu/data/"

    if len(sys.argv) > 2:
        data_folder = sys.argv[1]

    input_sample = open("/home/ubuntu/circrna/samples.csv", "w")
    input_sample.write("Sample_ID,Read1,Read2,Bam\n")
    for _, _, filenames in os.walk(data_folder):
        filenames.sort()
        for i, f_name in enumerate(filenames):
            line = "YY" + str(i) + "," + data_folder + f_name + ",NA,NA\n"
            input_sample.write(line)

    input_sample.close()
