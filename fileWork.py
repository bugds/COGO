import pickle
import codecs

def pickleFile(file):
    pickled = codecs.encode(pickle.dumps(file), 'base64').decode()

    return pickled

def unpickleFile(file):
    unpickled = pickle.loads(codecs.decode(file.encode(), 'base64'))

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

def makeRadioButtons(someDict):
    buttons = ''
    keyTemplate = '''
            <p><fieldset id="{}">{}</p>{}
            </fieldset>'''
    valueTemplate = '''
                <p><input type="radio" value="{}" name="{}"{}>{}</p>'''

    for key, d in someDict.items():
        comment = ''
        first = True
        for k, v in d.items():
            if first:
                comment += valueTemplate.format(k, key, ' checked', k + ' ' + v)
                first = False
            else:
                comment += valueTemplate.format(k, key, '', k + ' ' + v)

        buttons += keyTemplate.format(key, key, comment)

    return buttons