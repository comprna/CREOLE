# To change this license header, choose License Headers in Project Properties.
# To change this template file, choose Tools | Templates
# and open the template in the editor.


#def the_process_function():
#    n = 20
#    for i in range(n):
#        time.sleep(1)
#        sys.stdout.write('\r'+'loading...  process '+str(i)+'/'+str(n)+' '+ '{:.2f}'.format(i/n*100)+'%')
#        sys.stdout.flush()
#    sys.stdout.write('\r'+'loading... finished               \n')
#    
#    
from progress.bar import Bar
bar = Bar('Processing', max=20)
for i in range(20):
    # Do some work
    bar.next()
bar.finish()

import time, sys
def loading():
    print("Loading...")
    for i in range(0, 100):
        time.sleep(0.1)
        sys.stdout.write(u"\u001b[1000D" + str(i + 1) + "%")
        sys.stdout.flush()
    print()

loading()


#~/Documents/tools/diamond makedb --in sirvs_ref.fa -d reference
