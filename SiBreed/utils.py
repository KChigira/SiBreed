from datetime import datetime
import sys

def time_stamp():
    return '[program:{}]'.format(datetime.now().strftime('%Y-%m-%d %H:%M:%S'))

def prepare_cmd(cmd):
    #this function make the command arguments splitted by only space.
    return ' '.join(cmd.split())

def call_log(out_dir, name, cmd):
    print(time_stamp(), 
          '!!ERROR!! {}\n'.format(cmd), 
          flush=True)
    print('please check {}/log/{}.log\n'.format(out_dir, name))

def command(out_dir):
    #Output command info
    command = ' '.join(sys.argv)
    fn = '{}/command.txt'.format(out_dir)
    with open(fn, 'w') as f:
        f.write('{}\n'.format(command))
