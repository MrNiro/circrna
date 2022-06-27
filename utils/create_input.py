import sys
import os


if __name__ == '__main__':
    data_folder = "/home/ubuntu/data/"

    if len(sys.argv) > 2:
        data_folder = sys.argv[1]

    input_sample = open("/home/ubuntu/circrna/samples.csv", "w")
    input_sample.write("Sample_ID,Read1,Read2,Bam\n")
    for _, dirs, _ in os.walk(data_folder):
        for d in dirs.sort():
            line = d + ","
            cur_path = data_folder + d
            for dpath, _, filenames in os.walk(cur_path):
                for f_name in filenames.sort():
                    if "MD5" in f_name:
                        filenames.remove(f_name)

                if len(filenames) == 2:
                    line += (cur_path + filenames[0] + ",")
                    line += (cur_path + filenames[1] + ",NA\n")
                else:
                    Exception("Missing files in " + cur_path)
                break

            input_sample.write(line)

    input_sample.close()
