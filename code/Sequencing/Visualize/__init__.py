def optional_ax(original_function):
    def possibly_new_ax(*args, **kwargs):
        ax_given = kwargs.get('ax')
        if not ax_given:
            fig, ax = plt.subplots(figsize=(12, 8))
            kwargs['ax'] = ax
        
        result = original_function(*args, **kwargs)
        
        figure_file_name = kwargs.get('save_as')
        if figure_file_name:
            assert not ax_given, 'Can\'t give ax and figure_file_name'
            fig.savefig(figure_file_name, bbox_inches='tight')
            plt.close(fig)
        
        return result
    
    return possibly_new_ax
