# Optimal cancer treatment using model predictive control
The cancer is one of the biggest healthcare challenges that the modern society faces. Although our ability to deal with this medical condition have significantly increased in the recent decades, we are still unable to treat all the kinds of cancer and many tumors, such as brain tumor called 'glioblastoma multiformis' have very high mortality.

Current clinical treatment protocols of cancer involve different treatment modalities, such as radiation therapy, chemotherapy, immunotherapy, and surgical treatment. All the different types of treatment have positive and negative aspects and various side effects associated with them. For example, while the radiation therapy can kill the significant number of cancerous within the whole tumor, it can also damage the surrounding tissues. Chemotherapy can kill cancerous cells quite effectively, but it is also effective in killing of the cells accross the whole body due to its systemic effects.

In this blog, I would like to demonstrate how the methods of systems biology and control engineering can help us to optimize the cancer treatment through scheduling of the chemotherapy and immunotherapy, using the deterministic model of the cancer proliferation. I also hope to demonstrate how the (bio)engineering approach and computational tools can be applied to cancer. 

## The biological model
First of all, we need to think about how we will describe the model of the cancerous growth mathematically. 
