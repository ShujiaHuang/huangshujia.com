---
title: '多Github账号的使用问题'
date: 2014-08-10 01:00:00+0800
categories:
    - 编程技术
tags:
    - github


---

同一台电脑，或者同一个（大型机）服务器账号下要使用多个Github账号，该咋办？

如果是单用户(single-user)，很方便，默认拿id_rsa与你的github服务器的公钥对比；如果是多用户（multi-user）如user1,user2,那么就不能用在user2的身上了，这个时候就要配置一下了：

**1、新建user2的SSH Key**  

```bash
#新建SSH key(若有需求，注意不要误删其他的id_rsa key)：
$ cd ~/.ssh # 切换到~/.ssh
ssh-keygen -t rsa -C mywork@email.com  # 新建SSH key
# 设置名称为id_rsa_hshujia (设置别的名字，防止因名字冲突而误删其他的id_rsa)
Enter filein which to save the key (~/.ssh/id_rsa): id_rsa_hshujia
```


**2、新密钥添加到SSH agent中**  

因为默认只读取id_rsa，为了让SSH识别新的私钥，需将其添加到SSH agent中：

添加至SSH agent

```bash
ssh-add  ~/.ssh/id_rsa_hshujia
```

如果出现Could not open a connection to your authentication agent的错误，就试着用以下命令：

```bash
ssh-agent bash && ssh-add ~/.ssh/id_rsa_hshujia
```

**3、修改config文件**

在~/.ssh目录下找到config文件，vi打开进行编辑，如果没有就创建：

```bash
# 创建config
touch config  
```

然后按如下形式修改：

```bash
# 该文件用于配置私钥对应的服务器
# Default github user(first@mail.com)**
Host github.com
    HostName github.com
    User "user1"    
    IdentityFile ~/.ssh/id_rsa 

    # second user(second@mail.com)
    # 建一个github别名，新建的帐号使用这个别名做克隆和更新
Host hshujia.github.com     
    HostName github.com
    User "user2"    
    IdentityFile ~/.ssh/id_rsa_hshujia
```

其规则就是：从上至下读取config的内容，在每个Host下寻找对应的私钥。这里将GitHub SSH仓库地址中的***git@github.com***替换成新建的Host别名如：***hshujia.github.com***，那么原地址是：**git@github.com**:user/Mywork.git，替换后应该是：**git@hshujia.github.com**:user/Mywork.git.

以下是我在（大型机）服务器上config文件的具体内容：
![My config](https://static.fungenomics.com/images/2021/03/fg.post.2014-08-30-Fig1.jpg-blog.fungenomics-20210327223424453.com)

**4、打开新生成的~/.ssh/id_rsa_hshujia.pub文件，将里面的内容添加到GitHub后台，此处默认大家懂得如何添加，就不详述了。**

完成之后，在终端命令行中测试：

```
ssh -T git@hshujia.github.com  # 检验key是否已被成功添加到Github后台
Hi hshujia! You've successfully authenticated, but GitHub does not provide shell access. # 若能看到类似于这样的一句话就说明已经成功添加。
```

**5、若以前是Global配置，则取消Global配置。**

因为`git pull` or `git push` 的时候识别的是邮箱，多个github账号，就有多个邮箱，我们自然不能使用global的user.email了。

```bash
# 1.取消global 
git config --global --unset user.name
git config --global --unset user.email

# 2.在各个对应repo的目录下，设置项目repo自己的user.email和user name
git config  user.email "xxxx@xx.com"
git config  user.name "xxx"
```

**6、应用（例子）**

写个Hello world 测试一下。打开<github.com>登陆自己的账号，并新建一个repo，我这里命名为hello-world，创建完成之后，我在命令行中的具体操作如下:

```bash
git init
git add README.md
git config user.name hshujia
git config user.email hshujia@qq.com
git commit -m "first commit"
git remote add origin git@hshujia.github.com:hshujia/hello-world.git # 注意修改git@后面的名字 
git push -u origin master
```

这里唯一需要注意的就是多了两个`git config`用于告知该repo是属于哪一个user的；还有就是`git remote add origin`后面的host名字需要做相应的修改，具体的情况也已经在上面的代码中指明了。 

*参考资料*  

> 1. [git初体验（七）多账户的使用](http://www.cnblogs.com/BeginMan/p/3548139.html)
> 2. [Git的多账号如何处理](https://gist.github.com/suziewong/4378434)  

欢迎关注我的个人公众号：**helixminer（碱基矿工）**

![helixminer-QRCode](https://static.fungenomics.com/images/2021/03/helixminer-mid-red-20210327223413233-20210327223424681.png)