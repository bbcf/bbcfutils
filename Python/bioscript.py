
try:
    import bsPlugins
except ImportError:
    print "bsPlugins library must be in your PYTHONPATH in order to use this command-line tool"
    raise

import tw2.bs as twb
import tw2.dynforms as twd
import tw2.forms as twf


import sys
import os
import optparse

usage = "bioscript [OPTIONS] <plugin>"
plusage = lambda name: "bioscript {} [PLUGIN PARAMETERS]".format(name)
description = "Command line interface to use bsPlugins library"

opts = (
    ('-l', '--list', {'help': 'List available plugins', 'action': 'store_true', 'dest': 'list', 'default': False}),
    ('-i', '--info', {'help': 'Get the help from a plugin', 'action': 'store_true', 'dest': 'info', 'default': False}),
    ('-d', '--desc', {'help': 'Get the plugin description', 'action': 'store_true', 'dest': 'desc', 'default': False}),
    ('-a', '--args', {'help': 'Get the plugin arguments', 'action': 'store_true', 'dest': 'args', 'default': False}),
)

PLUGINS = list(bsPlugins.PLUGINS_FILES)
sorted(PLUGINS)

AVAILABLES = 'Available plugins are: {}{}{}'.format(os.linesep, os.linesep, ', '.join(PLUGINS))

DISPLAYS = {
    'description': lambda plugin: 'DESCRIPTION : {0}{1}'.format(os.linesep, plugin.info['description']),
    'inputs': lambda plugin: 'INPUTS : {0}{1}'.format(os.linesep, displ(plugin.info['in'], plugin=plugin)),
    'outputs': lambda plugin: 'OUTPUTS : {0}{1}'.format(os.linesep, displ(plugin.info['out']))
}


def displ(l, plugin=None):
    ids = ('id', 'type')
    with_options = ('radio', 'assembly', 'list')
    c = list(l)
    sorted(c, key=lambda d: d['id'])
    parameters = []
    for item in c:
        principal = '{item[id]} ({item[type]})'.format(item=item)
        second = ' '.join(['[{}]'.format(k, v) for k, v in item.iteritems() if k not in ids])
        third = ''
        if item['type'] in with_options:
            for att in plugin.info['output']().child.children:
                if att.id == item['id']:
                    third = att.options
                    third = ', '.join([x if isinstance(x, basestring) else x[0] for x in att.options])
        parameters.append('  {} {} {}'.format(principal, second, third))
    return os.linesep.join(parameters)


def main():
    parser = optparse.OptionParser(usage=usage, description=description)
    for opt in opts:
        parser.add_option(opt[0], opt[1], **opt[2])

    opt, args = parser.parse_args()

    # show the list of AVAILABLES plugins
    if opt.list:
        print AVAILABLES
        return 0

    # show help
    if len(args) == 0:
        parser.print_help()
        return 0

    # load the plugin
    plugin = load_plugin(args[0])
    if not plugin:
        return

    if opt.desc:
        print 'Usage: {1}{0}{0}{2}'.format(os.linesep, plusage(plugin.__module__.split('.')[-1]),
                                           DISPLAYS['description'](plugin))
        return 0
    elif opt.args:
        print 'Usage: {1}{0}{0}{2}{0}{0}{3}'.format(os.linesep, plusage(plugin.__module__.split('.')[-1]),
                                                    DISPLAYS['inputs'](plugin),
                                                    DISPLAYS['outputs'](plugin))
        return 0
    # show the help of the plugins specified
    elif opt.info or len(args) < 2:
        print 'Usage: {1}{0}{0}{2}{0}{0}{3}{0}{0}{4}'.format(os.linesep, plusage(plugin.__module__.split('.')[-1]),
                                                             DISPLAYS['description'](plugin),
                                                             DISPLAYS['inputs'](plugin),
                                                             DISPLAYS['outputs'](plugin))
        return 0
    else:
        d = parse_plugin_args(plugin, args[1:])
        p = plugin()
        try:
            result_plugin = p(**d)
            ret = "Plugin return {}".format(result_plugin)
            ofiles = "Output file(s) : {}".format(', '.join(['{} ({})'.format(path, t) for path, t in p.output_files]))
            print "{1}{0}{2}".format(os.linesep, ret, ofiles)
        except:
            print "Error in the plugin"
            raise
        finally:
            if len(p.tmp_files) > 0:
                print "Temporary directories(s) : {}".format(', '.join([t for t in p.tmp_files]))
                print
        return 0


def parse_plugin_args(plugin, args):
    accepted = [info['id'] for info in plugin.info['in']]
    # parse arguments list
    splitted_args = [a.split("=") for a in args]
    args = [(key, values.split(',')) if len(values.split(',')) > 1 else (key, values) for key, values in splitted_args]
    for arg_id, arg_args in args:
        if arg_id not in accepted:
            raise Exception("Argument '{}' is not accepted.".format(arg_id))

    # check if all requires args are presents
    required = [info['id'] for info in plugin.info['in'] if info.get('required', False)]
    for req in required:
        if req not in [i[0] for i in args]:
            raise Exception("Missing argument '{}'".format(req))

    # multiple fields must be arranged
    multiple = dict([(info['id'], info.get('multiple', False)) for info in plugin.info['in']])
    parsed = dict(args)
    arranged = dict(args)
    for k, v in parsed.iteritems():
        if multiple[k]:
            arranged[multiple[k]] = {k: v}
            del arranged[k]
    return arranged


def load_plugin(name):
    if name not in PLUGINS:
        print 'Plugin "{}" not in the plugin list. {}'.format(name, AVAILABLES)
        return
    import inspect
    import imp
    import traceback
    try:
        fp, pathname, description = imp.find_module(name, bsPlugins.__path__)
        try:
            p = imp.load_module(name, fp, pathname, description)
            clsmembers = inspect.getmembers(p, inspect.isclass)
            for name, clz in clsmembers:
                if clz.__module__ == p.__name__ and hasattr(clz, 'bs_plugin') and getattr(clz, 'bs_plugin') == 'bs-operation':
                    return clz
        except Exception as e:
            print 'Module "{}" not loaded cause : {}'.format(name, str(e))
            exc_type, exc_value, exc_traceback = sys.exc_info()
            print os.linesep.join(traceback.format_tb(exc_traceback))
        finally:
            fp.close()
    except ImportError as e:
        print 'Module {} not found.'.format(name)

if __name__ == '__main__':
    sys.exit(main())
