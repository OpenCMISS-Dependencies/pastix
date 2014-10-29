from cwrap.config import Config, File
import sys

print sys.argv
if len(sys.argv) > 1:
    files = [File(f) for f in sys.argv[2:]]
else:
    files = [File('test.h')]

if __name__ == '__main__':
    #config = Config('gccxml', files=files, save_dir = 'tests/result_gccxml')
    #config.generate()
    
    print '------------------------'
    print

    config_clang = Config('clang', files=files,
                          save_dir = sys.argv[1],
                          #include_dirs = [], #/usr/include/c++/4.2.1
                          #extern_name = '',
                          #implementation_name = '',
                          #language = 'c++',
                          )
    config_clang.generate()
