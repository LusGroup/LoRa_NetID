# LoRaPHY
这是一个在Matlab平台上实现LoRa解码、编码的代码。基于"[王继良老师团队](https://github.com/jkadbear/LoRaPHY)"的代码进行改进。

# 使用说明
+ LoRaPHY.m 文件包含了LoRa的调制、解调过程。并且实现了一些用于 debug 的画图、读取iq信号等基本功能。
+ general_lora_packet.m 文件实现了生成一个**特殊前导码**的 LoRa 数据包
+ test_sync_method.m 文件用于测试能否解码具有**特殊前导码**的 LoRa 数据包