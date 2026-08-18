# WASTER-Lite在AC自动机阶段使用的错配宽容搜索算法

实现方法：记录上一次错配到现在的位置。

例如：使用panama、nanako、banana组成的AC自动机宽容搜索banamako（每个word允许1次错配）：

null(inf) -1-> p(0) -> pa(1) -> pan(2) -> pana(3) -> panam(4) -> panama(5)<hit: panama> -> (fail) null -> end
          -2-> 
          -3->
          -4-> null(0) -> null