#!/usr/bin/env python

import sys
import os
import traceback
import argparse
import time
import pydevd
import smtplib
import datetime
import shutil
from jinja2 import Environment, FileSystemLoader
from essay import Essay

# gluster = os.environ['PROJECTS']
# gluster_2 = os.environ['PROJECTS2']
# gluster_done = os.environ['PROJECTS_DONE']
# gluster_2_done = os.environ['PROJECTS2_DONE']


# def project_exists(project_code):
#     gluster_project_path = os.path.join(gluster, project_code)
#     gluster_done_project_path = os.path.join(gluster_done, project_code)
#
#     gluster_2_project_path = os.path.join(gluster_2, project_code)
#     gluster_2_done_project_path = os.path.join(gluster_2_done, project_code)
#
#     if os.path.isdir(gluster_project_path):
#         return 1, gluster_project_path, gluster_done_project_path
#
#     elif os.path.isdir(gluster_2_project_path):
#         return 2, gluster_2_project_path, gluster_2_done_project_path
#
#     else:
#         sys.exit('ERROR: {0} project not found in Gluster/Gluster2.'.format(project_code))

# def convert_unicode_to_str(json_input):
#     if isinstance(json_input, dict):
#         return {convert_unicode_to_str(key): convert_unicode_to_str(value) for key, value in json_input.iteritems()}
#     elif isinstance(json_input, list):
#         return [convert_unicode_to_str(element) for element in json_input]
#     elif isinstance(json_input, unicode):
#         return json_input.encode('utf-8')
#     else:
#         return json_input
def get_mapping_stats(mapping_performance_file):
    with open(mapping_performance_file) as stats_file:
        next(stats_file)  # Skip header line
        lines = filter(None, (line.rstrip() for line in stats_file))  # Filter blank lines from stats file
        for line in lines:  # Iterate over non-blank lines in lines list ;)
            mapping_stats = line.split('\t')
            return mapping_stats[:-2]  # Return mapping_stats list, skip reads on target and % of reads on target.


def jinja_html_report(sample_list, mapping_stats_sample_list, seq_platform, sequencer, size, output_dir):
    #Get script path to load HTML template
    script_path = os.path.dirname(__file__)

    #Set report date time format
    now = datetime.datetime.now()
    date = now.strftime('%d-%m-%Y %H:%M')
    signature_date = now.strftime('%B %Y')

    #Jinja2 load path, template and variables
    env = Environment(loader=FileSystemLoader(script_path))
    template = env.get_template('report_template.html')
    res = template.render(sample_list=sample_list, mapping_stats_sample_list=mapping_stats_sample_list,
                          seq_platform=seq_platform, sequencer=sequencer, size=size,
                          date=date, signature_date=signature_date).encode('utf-8')

    #Output paths for CSS and pictures.
    css_output_path = os.path.join(output_dir, 'css')
    img_output_path = os.path.join(output_dir, 'img')

    #Copy css folder
    if os.path.exists(css_output_path):
        shutil.rmtree(css_output_path)
    css_src = os.path.join(script_path, 'css')
    shutil.copytree(css_src, css_output_path)

    #Copy img folder
    if os.path.exists(img_output_path):
        shutil.rmtree(img_output_path)
    img_src = os.path.join(script_path, 'img')
    shutil.copytree(img_src, img_output_path)

    #Open, write and close html
    html_output = output_dir + '/results_report.html'
    #pdf_output = '/share/gluster/tests/gmarco/informe_html' + '/results_report.pdf'
    f = open(html_output, 'w')
    f.write(res)
    f.close()

    print 'Generated HTML report: {0}'.format(html_output)


def print_logo():
    print '''   ___  _
  / _ \(_)__ ____ ___  ___
 / ___/ / _ `/ -_) _ \/ _ \\
/_/  /_/\_, /\__/\___/_//_/
        /__/
'''


# def check_results(project_dir):
#     pattern = '*.psv'
#
#     for root, dirs, files in os.walk(project_dir):
#         for filename in fnmatch.filter(files, pattern):
#             psv_path = (os.path.join(root, filename))
#             if os.stat(psv_path).st_size != 0:
#                 continue
#             else:
#                 print '{0} is empty, all the project analysis are not done yet !'.format(filename)


# def get_sample_list(project_dir):
#     sample_list = []
#     pattern = 'Sample_*'
#
#     rawdata_project_dir = os.path.join(project_dir, 'rawdata')
#     for root, dirs, files in os.walk(rawdata_project_dir):
#         for dirname in fnmatch.filter(dirs, pattern):
#             if dirname == 'Sample_Undetermined':
#                 continue
#             else:
#                 sample_list.append(re.search('Sample_(.*)', dirname).group(1))
#
#     return sample_list
def output_mgmt(output_path):
    """Creates output directories, even if nested in
    case that they don't exist."""

    if output_path == '.':
        raise ValueError(
            'The path {0} is not a valid option. Leave output option blank or specify a valid output file name.'.format(
                output_path))

    if not os.path.exists(output_path):
        print "Creating output directory..."
        os.makedirs(output_path)

    return os.path.abspath(output_path)


def send_mail(project_name, sample_list):
    sender = 'genetonic@sistemasgenomicos.com'
    #receivers = [ os.environ['USERMAIL'] ]
    to = [os.environ['USERMAIL']]
    #gmail_user = 'user@gmail.com'
    #gmail_pwd = 'pass'
    #smtpserver = smtplib.SMTP("smtp.gmail.com", 587)
    #smtpserver.ehlo()
    #smtpserver.starttls()
    #smtpserver.ehlo
    #smtpserver.login(gmail_user, gmail_pwd)
    #
    #header = 'To:' + to + '\n' + 'From: ' + gmail_user + '\n' + 'Subject:testing \n'
    #print header
    #msg = header + '\n this is test msg from genetonic \n\n'

    message = 'From: Genetonic Bioinformatic cluster <genetonic@sistemasgenomicos.com>\n'
    message += 'To: {0}\n'.format(to)
    message += 'Subject: Entrega de resultados: {0}\n'.format(project_name)
    message += 'Hola chic@s,\n\n' \
               'Ya se encuentran disponibles los resultados del provenientes del analisis {0} en el ' \
               'repositorio de datos.\n\n' \
               'Las muestras que se incluyen en este run son:\n\n'.format(project_name)

    for sample in sample_list:
        message += '- {0}\n'.format(sample)

    smtpobj = smtplib.SMTP('localhost')
    #smtpobj.sendmail(sender, to, message)
    print 'INFO: Successfully sent email'


def sample_stats_json_reader(essay_path):
    #Create essay object
    essay = Essay()

    #Read essay
    essay.read_essay(essay_path)

    #Get sample list
    essay_sample_id_list = sorted(essay.get_samples())

    #Create empty sample_structure
    mapping_stats_sample_list = []

    for sample in essay_sample_id_list:
        mapping_performance_file = essay.get_arg(info_type='tree', folder_type='analysis',
                                                 keys=['stats', '{0}'.format(sample), 'replicate1',
                                                       'mapping_performance.txt'])

        seq_platform = essay.get_arg(info_type='samples', folder_type='{0}'.format(sample),
                                     keys=['replicates', 'replicate1', 'seq_platform'])

        sequencer = essay.get_arg(info_type='samples', folder_type='{0}'.format(sample),
                                  keys=['replicates', 'replicate1', 'sequencer'])

        size = essay.get_arg(info_type='samples', folder_type='{0}'.format(sample),
                             keys=['replicates', 'replicate1', 'size'])

        mapping_stats_sample_list.append(
            {sample: {'mapping_performance': mapping_performance_file, 'mapping_stats': []}})

    for sample_item in mapping_stats_sample_list:
        for sample_id in sample_item.keys():
            mapping_performance_file = sample_item.get(sample_id).get('mapping_performance')

            mapping_stats = sample_item.get(sample_id).get('mapping_stats')
            mapping_stats.extend(get_mapping_stats(mapping_performance_file))

    return essay_sample_id_list, mapping_stats_sample_list, seq_platform, sequencer, size


def main():
    global args

    #Print ASCII logo
    print_logo()

    #Output dir management
    output_dir = output_mgmt(args.output_dir)

    #Read stats from Json
    essay_sample_id_list, mapping_stats_sample_list, seq_platform, sequencer, size = sample_stats_json_reader(
        args.essay_path)

    #We assume that all samples belong to same seq_plataform, sequencer and size.
    jinja_html_report(essay_sample_id_list, mapping_stats_sample_list, seq_platform, sequencer, size, output_dir)

    #send_mail(args.bf_project_name, sample_list)
    print 'DONE'


if __name__ == '__main__':
    try:
        start_time = time.time()

        parser = argparse.ArgumentParser(
            description='Creates HTML report for an specific exome essay.',
            version='alpha version')

        #parser.add_argument('-bf', dest='bf_project_name', action='store',
        #                    help='Bioinformatics project directory name',
        #                    required=True)

        parser.add_argument('-e', action='store', dest="essay_path",
                            help='Path to essay folder containing state.json file', required=True)

        parser.add_argument('-o', action='store', dest='output_dir',
                            help="Output directory.", required=True)

        parser.add_argument('-d', dest='debug', action='store_true',
                            help="Remote debug action.", required=False)

        args = parser.parse_args()
        if args.debug:
            user_ip = os.environ['USERIP']
            pydevd.settrace(user_ip, port=58484, stdoutToServer=True, stderrToServer=True)

        main()

    except SystemExit, e:  # sys.exit()
        #print 'System unexpected exit detected.'
        raise e

    except Exception, e:
        print 'ERROR: UNEXPECTED EXCEPTION'
        #print 'Removing DONE destination copy directory.'
        print str(e)
        traceback.print_exc()
        sys.exit(1)