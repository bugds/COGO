import pickle
import codecs

def pickleFile(obj):
    pickled = codecs.encode(pickle.dumps(obj), 'base64').decode()

    return pickled

def unpickleFile(obj):
    unpickled = pickle.loads(codecs.decode(obj.encode(), 'base64'))

    return unpickled

def makeOptionList(someDict):
    options = ''
    template = '''
                <option value="{}">{}</option>'''

    for key, value in someDict.items():
        if value == 1:
            comment = key + ': {} gene'.format(str(value))
        else:
            comment = key + ': {} genes'.format(str(value))

        options += template.format(key, comment)

    return options