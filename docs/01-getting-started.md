---
title: "01 Getting Started"
author: "Ines Lucia Patop: inespatop@brandeis.edu"
output: 
  html_document:
    toc: true
    toc_depth: 2
    number_sections: true
    df_print: paged
    theme: united
    toc_float: true
    collapsed: false
---



# Introduction to R and RStudio

## Objectives

After this section you should be able to:

1. Create projects and scripts in RStudio 
2. Install and load packages
3. Do basic statistical analysis and plots
4. Create basic for loops and if conditionals

## Introduction

We will be working with R and Rstudio, a user-friendly platform to use R. So, first you will need to download both [R](https://www.r-project.org/ ) and [RStudio](https://www.rstudio.com/). 

R is open-source and free programing language. This means that everybody, and I mean all of us, can collaborate, create packages and share things with others. This also means that everything we make in R is for public use and that, vice versa, we can use everything that anybody did in R. 

R is “object-orientated” a programing language. This basically means that most of the “things” in the environment are objects. I will not go deeper into it, but you can read more [here](https://cran.r-project.org/doc/contrib/Lemon-kickstart/kr_dobj.html). You will see in this class different types of objects. Each type of object has different properties which means we can not do everything with each type of object.

In this class we will give you some data and simple tasks so you can start playing around and get use to basic commands and operators. But first, we will need to establish some common language:

Let’s try to see some of the definitions more commonly used so we can understand them. I recommend you to do some research on your own.

* *Shell or Terminal*: the computer shell is a user interface to access to an operating system's services. In linux and mac is easier to access. You can read more [here](https://en.wikipedia.org/wiki/Shell_(computing)).

* *Programing language*: is a language in a wide sense, ie, it is a set of rules in which you can give orders to the computer. A computer basically, computes. Yes, that is why they are called computers. Then, you can use your computer as a calculator, as a table manager, to write text, etc. The idea is that if you know the language of a programing language you can tell your computer to do stuff for you. We will use R language. 

* *Script*: a text file that has the code (text in the programing language you are using) you will execute in the computer. In order for this to be executed you need to “copy-paste” the parts of code you wanna run into the shell. We will work in R-studio so this is done automatically with the ctr-enter command.

* *Working Directory*: in which folder (or directory) of your computer you will be running the code. As default R will run in your base-directory. This is the core of your computer. I recommend you to have a specific folder inside the documents folder with each of you projects and to run everything there. We will learn also how to set this up.

* *Environment*: All the objects, functions and everything you have loaded at that moment. This is a short-term memory thing. If we do not save this, everything we did will be deleted after we close the program. It would be like the words you have a in a text file if you do not save the file they will be lost.

* *Functions*: as mathematical functions, this takes inputs and generates outputs. We will see some in this class.

* *For loop, if, else*: this are basic logic operators that allows us to generate specific outputs. We can use it inside functions or as independent.

* *Package*:  R has the option to load many premade functions. The R package is a set of documented functions that someone did and put to be available for the rest of the community. We will use in this course many of them like: [DeSeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html), [ggplot2](https://ggplot2.tidyverse.org), etc.

## RStudio world

R Studio has everything integrated. You can have the script, the environment, the actual shell and a useful window to access data and to view the plots and see the manual and description of every function and package.

Here we can create a project in which everything will be integrated. This might not be the most efficient thing in terms of memory use but is the simplest for now. When you start to understand the logic of the programing language I encourage you to run R in different supports like typing “R” in your computer shell and run directly from command line.

Regarding RStudio see the image bellow:

<div class="figure" style="text-align: center">
<img src="./images/rstudio.png" alt="**CAPTION THIS FIGURE!!**" width="500px" />
<p class="caption">(\#fig:unnamed-chunk-1)**CAPTION THIS FIGURE!!**</p>
</div>

Upper left panel: the script. Everything you write here you can just execute by pressing CTR + ENTER. 

Upper right panel: the environment.

Down left panel: the “console”, work exactly as the terminal. Is the representation of the terminal just for R. If you want, you can run other things that are not R in the shell (not the console), if you are interested read this nice article https://support.rstudio.com/hc/en-us/articles/115010737148-Using-the-RStudio-Terminal 

Down Right: The multi panel that allows you to see and browse the files, plots and help pages.


## Create a project

1. Open RStudio:

Create a project: File > New  Project... > Existing  Directory > choose the folder. 

Automatically this should change the working directory to this folder. However, it is nice to check this we should run the following command:

```
getwd()  #Get working  directory
```

2. Create the script: File > New  File... >R Script  

Now we can start putting things in the script and running them in the “console” or representation of the shell in RStudio. So, let’s check the working directory and change it if needed.

To run the command, just write “getwd()” and press CTR+ENTER.

This should generate this in the “Console”:

```
getwd()
[1] "/Users/skl/Documents/Class/0_class" #here your working directory
```

Now we know how to run things in RStudio.

3. Save the script and the working directory: when you close the project it will automatically ask. 

I recommend to save at least the script everytime you can. File > save

## Useful shortcuts in RStudio:

* tab: auto-complete function

* control + the up arrow or command + up arrow: auto-complete tool that works only in the interactive console

* control + enter or command + enter: executes the selected lines of code

* control + s : save

## Installing and loading packages

When we open R the basic package with all the basic functions is loaded but we will use other functions from other packages. For now, we will start with ggplot2. This is not installed in our computer yet so we have to install it. This is done only ONCE in the computer. Then we need to load it in the current environment. This is done EVERY TIME we reopen R.

```
install.packages("ggplot2") #only once

library("ggplot2") #everytime we reopen R session
```

You might have noticed that I use “#” this a good way to add comments to the code. Any line that begins with a “#” will NOT be executed. 

We will be using packages that are part of [Bioconductor](https://bioconductor.org). This is a big repository for biological relevant packages. They are installed in this way:

```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager") #only once

BiocManager::install("DESeq2") #only once

library("DESeq2")

```

## Creating objects

Objects in R are created using the assign symbols: `<-` and `=` as follows: `object.name <- Object.Assingnemet`. 

This means that objects names should start with a letter and should NOT contain spaces. You can replace them with a dot or an underscore. 

Objects can be of different class and will be overwritten without any warnings.

If we execute the object name we will access it as better as the Console can do it. For numbers it is simple but you will notice that for other objects it is not.


```r
x<-0.5
x
class(x) #this tells us the class of the x object
```

```
## [1] 0.5
## [1] "numeric"
```

If we assign other thing to `x`, it will  be overwritten:


```r
x<-1 #this is overwriting x
x
class(x) #this tells us the class of the x object
```

```
## [1] 1
## [1] "numeric"
```

We can create as many objects as we want

```r
y<-2.3
y
class(y)
```

```
## [1] 2.3
## [1] "numeric"
```

## R as a calculator

R will operate as a calculator for numbers. It has a lot of prebuilt functions. Let’s see one obvious.


```r
x+y
(x+y)/2
```

```
## [1] 3.3
## [1] 1.65
```

We can now apply this sum a new object


```r
z<-x+y 
z
```

```
## [1] 3.3
```

And divide it by 2. 


```r
z/2
```

```
## [1] 1.65
```

Which is the same as doing the mean of x and y.


```r
mean(c(x,y))#what is this c()????
```

```
## [1] 1.65
```

To create a vector of numbers (or anything), you can just use c(n1,n2), you can also store this vector in a new object.


```r
v<-c(x,y)
v
class(v)
```

```
## [1] 1.0 2.3
## [1] "numeric"
```

You can now use it inside fucntions.

```r
mean(v)
```

```
## [1] 1.65
```

We can add (append) more elemts to the vector.


```r
t<-c(v,5)
t
mean(t)
#...and so on
```

```
## [1] 1.0 2.3 5.0
## [1] 2.766667
```
## Functions

We can think as functions exactly as we know mathematical functions. They take an input and generate an output. 

$$
y=f(x)
$$
$$
output=function(input)
$$
In R language that looks like this: `nameofunction<-function(x){ what the function does }` and the can do literally anything R can do. Math, plot, modify tables, etc.

Lets create a new mean function. We will call is mean.us. This will take the mean of two numbers. They will be called `n1` and `n2`.


```r
mean.us<-function(n1,n2){
  y<-((n1+n2)/2)
  return(y) #return is the one part of this function that actually makes it to return a value
}

#lets see if this works
mean.us(2,3)
```

```
## [1] 2.5
```
We can make it even more fancy and print a message


```r
mean.us<-function(n1,n2){
  y<-((n1+n2)/2)
  return(paste0("The mean of ",n1," and ",n2," is: ",y))
}

mean.us(2,3)
```

```
## [1] "The mean of 2 and 3 is: 2.5"
```
We can do now the mean plus 1

```r
mean.plus1<-function(n1,n2){
  y<-((n1+n2)/2+1)
  return(paste0("The mean of ",n1," and ",n2," plus one is: ",y))
}

mean.plus1(2,3)
```

```
## [1] "The mean of 2 and 3 plus one is: 3.5"
```

## Loops

Loops are useful to apply a function or an action to multiple objects. We will see `for loops` but be aware that another common loop type is the `while loop`.

*For loops* will go over each element of a vector or list provided. In R language they look like this: `for (x in vector) { DO SOMETHING }`. And it literally means that it will go over each element of the vector, each time x will take the value of the element it goes that time and will do something. 


```r
#i this is a little complex so lets go one step at a time: i is an object that will be getting the value of each elemnt we go thru

for(x in c(1,2,3,4)){
  print(paste0("x is: ",x))
}
```

```
## [1] "x is: 1"
## [1] "x is: 2"
## [1] "x is: 3"
## [1] "x is: 4"
```
We can now do something more complex inside the loop.


```r
for(x in c(1,2,3,4)){
  c=x/2
  print(paste0("c is: ",c))
}
```

```
## [1] "c is: 0.5"
## [1] "c is: 1"
## [1] "c is: 1.5"
## [1] "c is: 2"
```

We can also loop over vectors that are already in our list of objects.


```r
#t is mande of many elements already
t
#if we want to sum 1 to each element in t and print it out we can do as follows
for(i in t){
  print(i)
  c=i+1
  print(c)
}
```

```
## [1] 1.0 2.3 5.0
## [1] 1
## [1] 2
## [1] 2.3
## [1] 3.3
## [1] 5
## [1] 6
```

## Conditions

R assigns using `=` and compares using `==`,`<` and `>`. It can then use `if` and `else` to generate conditions and actions. This can be more complex by adding AND, OR gates with `&` and `|` respectively.


```r
x=1
#lets explore x
x
```

```
## [1] 1
```
Compare, equal to:

```r
x==1
```

```
## [1] TRUE
```
Smaller than:

```r
x<2
```

```
## [1] TRUE
```
Bigger than:

```r
x>2
```

```
## [1] FALSE
```

Apply this comparisons to an if/else:


```r
if (x < 2){
  c=x+1
  print(paste0("c is: ", c))
} else {
  print("x is too big")
}
```

```
## [1] "c is: 2"
```
Change the condition:
 

```r
if (x < 1){
  c=x+1
  print(paste0("c is: ", c))
} else {
  print("x is too big")
}
```

```
## [1] "x is too big"
```

## Character objects

So far, we saw opperations with numbers. But R can have objects with words or characters.

We can have a complex phrase

```r
phrase<-"I am a beatufill phrase. Hello world"
phrase
```

```
## [1] "I am a beatufill phrase. Hello world"
```
Or a vector of character elements:

```r
chrvec<-c(phrase,"abcde","67")
chrvec
```

```
## [1] "I am a beatufill phrase. Hello world"
## [2] "abcde"                               
## [3] "67"
```


```r
class(chrvec)
```

```
## [1] "character"
```
## Identify

We can operate over them. For example, we can select the elemnt that has certain characteristic. For this we will use [grep](https://www.rdocumentation.org/packages/base/versions/3.6.2/topics/grep). This is based on [regualar expression](https://rstudio.com/wp-content/uploads/2016/09/RegExCheatsheet.pdf).

Lets select the elemnts that have the letter d.


```r
grep(pattern = "d",x = chrvec,value = T)
```

```
## [1] "I am a beatufill phrase. Hello world"
## [2] "abcde"
```


```r
grep(pattern = "d",x = chrvec,value = F) #What changed between the previous one and this one?
```

```
## [1] 1 2
```
Clearly, the fist two elements have the letter d. What if we want to select the one that has it at the end?

```r
grep(pattern = "d$",x = chrvec,value = T)
```

```
## [1] "I am a beatufill phrase. Hello world"
```

The same can be done for the begining. This of course gives us an empty result.


```r
grep(pattern = "^d",x = chrvec,value = T)
```

```
## character(0)
```
## Modify

We can modify the objects. We will use [gsub](https://stackoverflow.com/questions/35655485/replace-with-space-using-gsub-in-r)


```r
gsub(pattern = "d",replacement = "I am replacing d",x = chrvec)
```

```
## [1] "I am a beatufill phrase. Hello worlI am replacing d"
## [2] "abcI am replacing de"                               
## [3] "67"
```

## Useful resources

* R for data science: https://r4ds.had.co.nz
* A package to learn R 

## Activity: 

**Create a new funciton that has inside a for loop and other that has an if/else condition.**



```r
plus.one.onlyifpos <- function(n){
  if(n > 0){
    return(n+1)
  } else {
    return("number is negative")
  }
}

plus.one.onlyifpos(20)
plus.one.onlyifpos(-20)
```
