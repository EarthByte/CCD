#!/usr/bin/env python3

# Run `python3 timescale_conversion.py -h` to find out how to use this script
# Checkout readme.md at https://github.sydney.edu.au/EarthByte/EarthByteWorkflows/blob/ff2c9b548732a95a4b0301b5daf9426cdf48f7de/TimescaleConversion/readme.md

import sys, csv, math, argparse
import xml.etree.ElementTree as ET

PRECISION = 0.05
MAX_AGE = 180
FLAGS = []

flag_desc_map = {
    0: "(no change) invalid line",
    1: "(no change) age exceeds user-defined max age",
    2: "(no change) absolute age",
    3: "(changed) not found in the current timescale, interpolated",
    4: "(changed) found in current timescale and the time in the new timescale is different",
    5: "(no change) the time in new timescale is the same with current time",
    6: "(no change) age exceeds the max age in the new timescale",
    7: "(no change) not found in the current timescale and the new interpolated time is the same with the current time",
    8: "(no change) not found in the current timescale, age exceeds the max age in the current timescale",
}

# add new timescale here
# you also need to update timescales.txt accordingly
timescale_names = [
    "Gee_Kent_2007",
    "Cande_Kent_1995",
    "Gradstein_2005",
    "Ogg_2012",
    "Ogg_2020",
]


def isfloat(x):
    """check if the input x is a floating-point number"""
    try:
        a = float(x)
    except ValueError:
        return False
    else:
        return True


def read_timescales(fn):
    """read in the timescales from timescales file, such as timescales.txt

    :param fn: the file name which contains the timescale data

    """
    timescales = []
    for _ in range(len(timescale_names)):
        timescales.append([])

    with open(fn, "r", newline="") as csvfile:
        reader = csv.reader(csvfile, delimiter="\t")
        for row in reader:
            if row[0].startswith("#"):
                continue
            if len(row) != len(timescale_names) * 4 + 2:
                print(f"invalid row: {row}")
                continue
            for idx in range(len(timescales)):
                try:
                    timescales[idx].append(float(row[idx * 2 + 1]))
                except ValueError:
                    timescales[idx].append(float("nan"))

                try:
                    timescales[idx].append(float(row[idx * 2 + 2]))
                except ValueError:
                    timescales[idx].append(float("nan"))

    with open("timescales.log", "w+", newline="") as ts_log:
        for idx in range(len(timescale_names)):
            ts_log.write(f"{timescale_names[idx]: <20}")
        ts_log.write("\n")
        if not len(timescales) > 0:
            print(f"ERROR: empty timescales {timescales}")
            sys.exit(-1)
        for idx_1 in range(len(timescales[0])):
            for idx_2 in range(len(timescale_names)):
                ts_log.write(f"{timescales[idx_2][idx_1]: <20}")
            ts_log.write("\n")

    return timescales


def get_max_age(ts):
    """get the max age from a timescale

    :param ts: a list of ages in a timescale

    """
    ret = 0.0
    for row in ts:
        if isfloat(row):
            ret = row
    return ret


def get_age_index(age, ts):
    """find the index value of an age from a timescale
    if an exact match is found, return (True, index)
    if no exact match is found, return (Fasle, index), the index is the first value greater than the input age

    :param age: the input age to be searched against
    :param ts: a list ages from a timescale

    """
    for i in range(len(ts)):
        if math.isnan(ts[i]):
            continue

        if abs(ts[i] - age) > PRECISION:  # not close enough
            if age > ts[i]:
                continue
            else:
                # print('age not found in the timescales.')
                return False, i
        else:  # found a match
            return True, i
    return False, len(ts)


def interp_age(age, start_idx, end_idx, old_ts, new_ts):
    """interpolate the new timescale to get a new age

    :param age: the input age
    :param start_idx: the start index of a timescale segment
    :param end_idx: the end index of a timescale segment
    :param old_ts: a list of ages from the old timescale
    :param new_ts: a list of ages from the new timescale
    """
    # find the non-nan boundaries
    while start_idx >= 0:
        if math.isnan(old_ts[start_idx]) or math.isnan(new_ts[start_idx]):
            start_idx -= 1
        else:
            break

    while end_idx < len(old_ts):
        if math.isnan(old_ts[end_idx]) or math.isnan(new_ts[end_idx]):
            end_idx += 1
        else:
            break

    if end_idx == len(old_ts):  # the age is too large for the new time scales
        return False, age
    # compute the new age
    ratio = (age - old_ts[start_idx]) / (old_ts[end_idx] - old_ts[start_idx])
    return True, new_ts[start_idx] + ratio * (new_ts[end_idx] - new_ts[start_idx])


def get_new_age(age, old_ts, new_ts):
    """get a new age from the new timescale
    if a new different age is computed, return (True, new_age)
    if no new age is computed, turn (False, old_age).
    This happens when the age in the new timescale is the same as the age in the old timescale.
    Or when the input age exceeds the upper limit in the new timescale.

    :param age: the input age
    :param old_ts: a list of ages from the old timescale
    :param new_ts: a list of ages from the new timescale

    """
    _, idx = get_age_index(
        age, old_ts
    )  # find out the location of the age in the old timescale
    if idx >= len(old_ts):  # out of boundary
        return False, age
    if abs(old_ts[idx] - age) < PRECISION:  # find a exact match
        if not math.isnan(new_ts[idx]):
            return False, new_ts[idx]
    assert idx >= 1
    return interp_age(age, idx - 1, idx, old_ts, new_ts)


def process_rotation_file(fn, new_fn, old_ts, new_ts):
    """replace the timescale with a new one in a rotation file

    :param fn: the input file name
    :param new_fn: the output file name
    :param old_ts: a list of ages from the old timescale
    :param new_ts: a list of ages from the new timescale

    """
    old_lines = []
    new_data = []
    not_found_age = []
    count = 0
    with open(fn, "r", newline="") as rot_file:
        for line in rot_file:
            count += 1
            line = line.lstrip()
            six_numbers = []
            comments = ""
            tmp = line.split("!", 1)  # separate the comments part
            if len(tmp) == 2:
                comments = tmp[1]

            six_numbers = tmp[0].split()
            if (
                len(six_numbers) < 6
                or not six_numbers[0].isdigit()
                or not isfloat(six_numbers[1])
            ):
                new_data.append(line)
                old_lines.append(line)
                FLAGS.append(0)  # no change, invalid line
                continue

            # save the old data
            new_line = ""
            for number in six_numbers:
                new_line += "{0: <10}".format(number)
            old_lines.append(new_line + " ! " + comments)

            # if it is absolute age, we don't change it
            # if age is greater than MAX_AGE, we do nothing
            age = float(six_numbers[1])
            # print(age)
            if comments.find("@absage") == -1 and age < MAX_AGE:
                found, idx = get_age_index(age, old_ts)
                if not found:
                    not_found_age.append(
                        "Line {0}:\n{1}\n".format(count, line)
                    )  # record the not-found-age

                # try the best to get a new age
                interp_flag, new_age = get_new_age(age, old_ts, new_ts)
                six_numbers[1] = round(new_age, 3)
                if str(six_numbers[1]) == str(age):
                    if get_max_age(old_ts) < age:
                        FLAGS.append(
                            8
                        )  # no change, age exceeds the max age in the current timescale
                    elif get_max_age(new_ts) < age:
                        FLAGS.append(
                            6
                        )  # no change, age exceeds the max age in the new timescale
                    else:
                        if interp_flag:
                            FLAGS.append(
                                7
                            )  # no change, the new interpolated time is the same with current time
                        else:
                            FLAGS.append(
                                5
                            )  # no change, the time in new timescale is the same with current time
                else:
                    if found:
                        FLAGS.append(4)  # changed, found in current timescale
                    else:
                        FLAGS.append(3)  # changed, interpolated

            else:
                if age >= MAX_AGE:
                    FLAGS.append(1)  # no change, exceed user-defined max age
                elif comments.find("@absage") != -1:
                    FLAGS.append(2)  # no change, absolute age

            new_line = ""
            for number in six_numbers:
                new_line += "{0: <10}".format(number)
            new_data.append(new_line + " ! " + comments)

    with open(new_fn, "w+", newline="") as new_rot_file:
        for row in new_data:
            new_rot_file.write(row)

    with open(
        "time_not_found_in_timescale.log", "w+", newline=""
    ) as not_found_ages_file:
        not_found_ages_file.write(
            "#The lines below contains the time values which are not found in the current timescale and "
            + "there is no @absage flag in the lines.\n\n"
            ""
        )
        for row in not_found_age:
            not_found_ages_file.write(row)

    with open("debug.log", "w+", newline="") as debug_file:
        for row in zip(old_lines, new_data, FLAGS):
            a_1 = row[0].split()[1]
            a_2 = row[1].split()[1]

            debug_file.write(flag_desc_map[row[2]])
            if a_1 != a_2:
                debug_file.write(
                    ". The age has been changed from {0} to {1}!\n".format(a_1, a_2)
                )
                if float(a_1).is_integer():
                    debug_file.write(
                        "(Warning) A whole number has been converted to the new timescale. Are you sure this is what you intended?\n"
                    )
            else:
                debug_file.write("\n")
            debug_file.write("old: " + row[0])
            debug_file.write("new: " + row[1])
            debug_file.write("\n")


def process_gpml_file(fn, new_fn, old_ts, new_ts):
    """replace the timescale with a new one in a gpml file

    :param fn: the input file name
    :param new_fn: the output file name
    :param old_ts: a list of ages from the old timescale
    :param new_ts: a list of ages from the new timescale

    """
    debug_info = []
    namespaces = {
        "gpml": "http://www.gplates.org/gplates",
        "gml": "http://www.opengis.net/gml",
        "xsi": "http://www.w3.org/XMLSchema-instance",
    }
    for ns in namespaces:
        ET.register_namespace(ns, namespaces[ns])
    # print(namespaces)
    tree = ET.parse(fn)
    root = tree.getroot()
    for member in root.findall("./gml:featureMember", namespaces):
        fid = member.find(".//gpml:identity", namespaces)
        fid = fid.text
        time = member.find(".//gml:begin//gml:timePosition", namespaces)

        if isfloat(time.text):
            age = float(time.text)
            if age > MAX_AGE or age > get_max_age(old_ts) or age > get_max_age(new_ts):
                continue

            # try the best to get a new age
            _, new_age = get_new_age(age, old_ts, new_ts)
            new_age = round(new_age, 3)
            if str(age) != str(new_age):
                debug_info.append(
                    "The age has been changed from {0} to {1} in feature ({2}).".format(
                        age, new_age, fid
                    )
                )
                time.text = str(new_age)

    tree.write(
        new_fn, short_empty_elements=False, encoding="utf-8", xml_declaration=True
    )
    with open("debug.log", "w+", newline="") as debug_file:
        for row in debug_info:
            debug_file.write(row + "\n")


def process_txt_file(
    input_filename,
    output_filename,
    old_timescale,
    new_timescale,
):
    """process text file. the input example is Site1410A_age_depth.txt.

    :param input_filename: the input file name
    :param output_filename: the output file name
    :param old_timescale: a list of ages from the old timescale
    :param new_timescale: a list of ages from the new timescale

    """
    output_data = []
    debug_info = []
    count = 0
    with open(input_filename, "r", newline="") as txt_file:
        for line in txt_file:
            count += 1
            if line.startswith(">"):  # comment line
                output_data.append(line)
                continue
            values = line.split()
            if len(values) < 1:  # empty line
                output_data.append(line)
                continue

            try:
                age = float(values[0])
                new_flag, new_age = get_new_age(age, old_timescale, new_timescale)
                if not new_flag:  # the new age is the same as the old age
                    output_data.append(line)
                else:
                    debug_info.append(
                        f"Line {count}: age changed from {age} to {round(new_age,6)}\n"
                    )
                    output_data.append(
                        f"{round(new_age,6):<15}{''.join(values[1:])}\n"
                    )  # write out the new age
            except ValueError:
                output_data.append(line)  # invalid age column
                continue

    with open(output_filename, "w+", newline="") as output_file:
        for row in output_data:
            output_file.write(row)
    with open("debug.log", "w+", newline="") as debug_file:
        for row in debug_info:
            debug_file.write(row)


if __name__ == "__main__":
    __description__ = """Convert the time values in a rotation file from one timescale to another.
    Timescales:
        Gee_Kent_2007 
        Cande_Kent_1995 
        Gradstein_2005
        Ogg_2012
        Ogg_2020
        
    The timescales.txt needs to be in the current working directory. Otherwise, you need to use -t to specify the timescale file.
        
    Example:
        python3 timescale_conversion.py -m 140 Global_EarthByte_230-0Ma_GK07_AREPS.rot output.rot Gee_Kent_2007 Gradstein_2005
        
        python3 timescale_conversion.py -m 140 Global_EarthByte_230-0Ma_GK07_AREPS_IsoCOB.gpml output.gpml Gee_Kent_2007 Gradstein_2005
        
        python3 timescale_conversion.py -m 140 -p 0.00001 Site1410A_age_depth.txt  output.txt Gee_Kent_2007 Gradstein_2005
        
    For rotation file, the script will create four(4) files. 
    1. the new rotation file
    2. time_not_found_in_timescale.log (ages not found in the current timescale and no @absage flag)
    3. timescales.log (timescales)
    4. debug.log (information about what has been done to each line)
    
    For GMPL file, the script will create three(3) files.
    1. The new gpml file with the new ages
    2. timescales.log (timescales)
    3. debug.log (information about which features' ages have been changed)
    
    """
    # The command-line parser.
    parser = argparse.ArgumentParser(
        description=__description__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        "-m",
        "--max-time",
        type=float,
        required=False,
        metavar="max_time",
        help="The time value which is greater than this max_time will not be processed.",
    )

    parser.add_argument(
        "-p",
        "--precision",
        type=float,
        required=False,
        metavar="precision",
        help="the precision tolerance which is used to decide if two float numbers are equal",
    )

    parser.add_argument(
        "-t",
        "--timescale",
        type=str,
        required=False,
        metavar="timescale",
        help="the timescale filename",
    )

    parser.add_argument(
        "input_filename",
        type=str,
        help="the name of the input file(can be rotation file or gpml file)",
    )

    parser.add_argument("output_filename", type=str, help="the name of the output file")

    parser.add_argument(
        "current_timescale",
        type=str,
        help="the name of the timescale which is currently being used in the rotation file",
    )

    parser.add_argument(
        "new_timescale",
        type=str,
        help="the name of the new timescale to which the time values will be converted",
    )

    # Parse command-line options.
    args = parser.parse_args()

    timescale_fn = "timescales.txt"

    if args.precision:
        PRECISION = args.precision
    if args.max_time:
        MAX_AGE = args.max_time
    if args.timescale:
        timescale_fn = args.timescale

    # read in the timescales
    timescales = read_timescales(timescale_fn)

    if args.current_timescale not in timescale_names:
        print(
            "Unknow timescale: {}. Do nothing and quit.".format(args.current_timescale)
        )
        sys.exit(1)
    if args.new_timescale not in timescale_names:
        print("Unknow timescale: {}. Do nothing and quit.".format(args.new_timescale))
        sys.exit(1)

    current_timescale_index = timescale_names.index(args.current_timescale)
    new_timescale_index = timescale_names.index(args.new_timescale)

    if args.input_filename.endswith(".rot"):
        process_rotation_file(
            args.input_filename,
            args.output_filename,
            timescales[current_timescale_index],
            timescales[new_timescale_index],
        )

        print("Done! \nThe new rotation file is {}.".format(args.output_filename))
        print(
            "The lines with ages not found in {} and no @absage flag are saved in time_not_found_in_timescale.log.".format(
                args.current_timescale
            )
        )
        print("The timescales are saved in timescales.log.")
        print("The full debug information is saved in debug.log.")
        print("Any problem, please contact michael.chin@sydney.edu.au")
    elif args.input_filename.endswith(".gpml"):
        process_gpml_file(
            args.input_filename,
            args.output_filename,
            timescales[current_timescale_index],
            timescales[new_timescale_index],
        )

        print("Done! \nThe new gpml file is {}.".format(args.output_filename))
        print("The debug information has been recorded in debug.log.")
    elif args.input_filename.endswith(".txt"):
        process_txt_file(
            args.input_filename,
            args.output_filename,
            timescales[current_timescale_index],
            timescales[new_timescale_index],
        )

        print("Done! \nThe new text file is {}.".format(args.output_filename))
        print("The debug information has been recorded in debug.log.")
    else:
        print("Unsupported input file: {}.".format(args.input_filename))
        sys.exit(1)
