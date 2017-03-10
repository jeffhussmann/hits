import os
import bokeh.models

def load_file(key, ext):
    fn = os.path.join(os.path.dirname(__file__), '{0}.{1}'.format(key, ext))
    with open(fn) as fh:
        callback = fh.read()
    return callback

def build_callback(key, js=False, format_kwargs=None, args=None):
    if args is None:
        args = {}
    if format_kwargs is None:
        format_kwargs = {}

    if js:
        ext = 'js'
        model = bokeh.models.CustomJS
    else:
        ext = 'coffee'
        model = bokeh.models.CustomJS.from_coffeescript

    code_template = load_file(key, ext)
    code = code_template.format(**format_kwargs)
    callback = model(code=code, args=args)

    return callback
