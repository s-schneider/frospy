import io
import pickle


class RenameUnpickler(pickle.Unpickler):
    def find_class(self, module, name):
        renamed_module = module
        if module == "frospy.setup.settings":
            renamed_module = "frospy.core.setup.settings"

        return super(RenameUnpickler, self).find_class(renamed_module, name)


def renamed_load(file_obj):
    return RenameUnpickler(file_obj).load()


def renamed_loads(pickled_bytes):
    file_obj = io.BytesIO(pickled_bytes)
    return renamed_load(file_obj)

x = '/quanta1/home/simons/splitting/modes/00t15-00s14/cross_coupling/TZ-Comp/deg20-20/setup.pickle'
with open(x, 'rb') as f:
    renamed_load(f)
