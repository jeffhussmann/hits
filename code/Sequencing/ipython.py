import ipywidgets

toggle = '''
<script>
code_show=true; 
function code_toggle() {
 if (code_show){
 $('div.input').hide();
 } else {
 $('div.input').show();
 }
 code_show = !code_show
} 
$( document ).ready(code_toggle);
</script>
To toggle code, click <a href="javascript:code_toggle()">here</a>.
'''

def binary_widget(default=True):
    values = [int(default), int(not default)]
    labels = [str(default), str(not default)]
    return ipywidgets.RadioWidget(values, labels=labels)
