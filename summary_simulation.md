This article provides analysis of the issue in the posterior predictive checking we discussed in <span class="timestamp-wrapper"><span class="timestamp">[2021-03-01 Mon] </span></span> meeting.


# TL;DR

It is advised to switch signs of the latent distance terms. A new model gives acceptable results for posterior predictive checking of the response type. The computational method I implemented doesn&rsquo;t seem to be wrong. So the issue might be coming from the statistical model itself.


# Models

The below is a list of models in consideration. They differ by the sign of the latent distance term. A baseline model is one with no latent distance term.

(-,+): (correct,incorrect)

\begin{align*}
h_{i1}(t) &= \lambda_{i1}(t) \exp(\beta_{i1} + \theta_{k1} - ||z_{k} - w_{i}||) \\
h_{i0}(t) &= \lambda_{i0}(t) \exp(\beta_{i0} + \theta_{k0} + ||z_{k} - w_{i}||) \\
\end{align*}

(+,-): (correct,incorrect)

\begin{align*}
h_{i1}(t) &= \lambda_{i1}(t) \exp(\beta_{i1} + \theta_{k1} + ||z_{k} - w_{i}||) \\
h_{i0}(t) &= \lambda_{i0}(t) \exp(\beta_{i0} + \theta_{k0} - ||z_{k} - w_{i}||) \\
\end{align*}


# Histogram of survival times

You need to open this document in Github using the web browser to see PDF files. The simulated times are divided into two panes based on the observed response type. The (-,+) model didn&rsquo;t have the issue of simulating the max survival time (31).

-   [(-,+) model](./chessB_pn/hist_time.pdf)
-   [(+,-) model](./chessB_pn/hist_time.pdf)
-   [baseline model](./chessB_no_latent/hist_time.pdf)


# Log Loss

Given the survival time and MCMC samples, the conditional probability of &ldquo;correct&rdquo; response is estimated as
\[
\hat p_{ik}^{(l)}(t) = \frac{h_{i1}^{(l)}(t)}{h_{i1}^{(l)}(t) + h_{i0}^{(l)}(t)},
\]
where \((l)\) denote the index of MCMC sample used to calculate the quantity.
We use \(Y_{ik} =
1\) to denote a correct response, and \(Y_{ik} = 0\) an incorrect response. Then, \(p_{i}^{(l)}(t)\) can be served as the prediction probability of \(Y_{ik}\), which allows using classification performance metrics. For easiness of evaluation, we use the log loss (over respondent and item)
\[
\sum_{l=1} ^{L}\left\{ y_{ik} \log \hat p_{ik}^{(l)}(t) + (1-y_{ik}) \log (1 - \hat p_{ik}^{(l)})\right\}.
\]

The log loss are concatenated over respondents for each item. Then, we calculated differences of concatenated log loss: 1) baseline - (-,+) log loss and 2) baseline - (+,-) log loss. So, the larger is the better. We provide a series of box plots for the differences over all items.

On average, the (-,+) model gives smaller log losses than the baseline for all correct responses, but and the (+,-) model fails to do so for both response types. In terms of &ldquo;posterior density&rdquo; and &ldquo;likelihood&rdquo;, these two models give similar fit, and also they show better fit than the baseline model.
The simulated survival time is used in the 1st PDF, and the observed survival time is used in the 2nd PDF. They show similar patterns (the results is a bit more promising with the simulation though), so the simulation is not likely to be a &ldquo;culprit&rdquo;. I also didn&rsquo;t find any error in the MCMC sampler. Perhaps, this could be what we get from switching signs of the latent distance.

-   [With simulated survival time](sim_logLoss_diff_by.pdf)
-   [With observed survival time](LogLoss_diff_by.pdf)

Furthermore, tables for means and standard errors for the differences (with simulated survival time) are listed as the following:


## Log loss differences

<table border="2" cellspacing="0" cellpadding="6" rules="groups" frame="hsides">


<colgroup>
<col  class="org-right" />

<col  class="org-right" />

<col  class="org-right" />

<col  class="org-right" />

<col  class="org-right" />
</colgroup>
<thead>
<tr>
<th scope="col" class="org-right">item</th>
<th scope="col" class="org-right">pos\<sub>neg</sub></th>
<th scope="col" class="org-right">se\<sub>pn</sub></th>
<th scope="col" class="org-right">neg\<sub>pos</sub></th>
<th scope="col" class="org-right">se\<sub>np</sub></th>
</tr align="right">
</thead>

<tbody>
<tr>
<td class="org-right">1</td>
<td class="org-right">-0.935</td>
<td class="org-right">0.039</td>
<td class="org-right">0.692</td>
<td class="org-right">0.041</td>
</tr align="right">


<tr>
<td class="org-right">2</td>
<td class="org-right">-0.944</td>
<td class="org-right">0.038</td>
<td class="org-right">0.636</td>
<td class="org-right">0.04</td>
</tr align="right">


<tr>
<td class="org-right">3</td>
<td class="org-right">-0.892</td>
<td class="org-right">0.04</td>
<td class="org-right">0.46</td>
<td class="org-right">0.04</td>
</tr align="right">


<tr>
<td class="org-right">4</td>
<td class="org-right">-0.938</td>
<td class="org-right">0.038</td>
<td class="org-right">0.511</td>
<td class="org-right">0.04</td>
</tr align="right">


<tr>
<td class="org-right">5</td>
<td class="org-right">-0.919</td>
<td class="org-right">0.039</td>
<td class="org-right">0.54</td>
<td class="org-right">0.04</td>
</tr align="right">


<tr>
<td class="org-right">6</td>
<td class="org-right">-0.882</td>
<td class="org-right">0.04</td>
<td class="org-right">0.326</td>
<td class="org-right">0.038</td>
</tr align="right">


<tr>
<td class="org-right">7</td>
<td class="org-right">-0.888</td>
<td class="org-right">0.039</td>
<td class="org-right">0.331</td>
<td class="org-right">0.037</td>
</tr align="right">


<tr>
<td class="org-right">8</td>
<td class="org-right">-0.854</td>
<td class="org-right">0.041</td>
<td class="org-right">0.221</td>
<td class="org-right">0.035</td>
</tr align="right">


<tr>
<td class="org-right">9</td>
<td class="org-right">-0.613</td>
<td class="org-right">0.045</td>
<td class="org-right">0.059</td>
<td class="org-right">0.029</td>
</tr align="right">


<tr>
<td class="org-right">10</td>
<td class="org-right">-0.771</td>
<td class="org-right">0.043</td>
<td class="org-right">0.35</td>
<td class="org-right">0.04</td>
</tr align="right">


<tr>
<td class="org-right">11</td>
<td class="org-right">-0.747</td>
<td class="org-right">0.043</td>
<td class="org-right">0.089</td>
<td class="org-right">0.031</td>
</tr align="right">


<tr>
<td class="org-right">12</td>
<td class="org-right">-0.606</td>
<td class="org-right">0.045</td>
<td class="org-right">0.034</td>
<td class="org-right">0.026</td>
</tr align="right">


<tr>
<td class="org-right">13</td>
<td class="org-right">-0.659</td>
<td class="org-right">0.043</td>
<td class="org-right">0.165</td>
<td class="org-right">0.033</td>
</tr align="right">


<tr>
<td class="org-right">14</td>
<td class="org-right">-0.593</td>
<td class="org-right">0.045</td>
<td class="org-right">0.049</td>
<td class="org-right">0.028</td>
</tr align="right">


<tr>
<td class="org-right">15</td>
<td class="org-right">-0.814</td>
<td class="org-right">0.042</td>
<td class="org-right">0.178</td>
<td class="org-right">0.033</td>
</tr align="right">


<tr>
<td class="org-right">16</td>
<td class="org-right">-0.328</td>
<td class="org-right">0.043</td>
<td class="org-right">-0.013</td>
<td class="org-right">0.024</td>
</tr align="right">


<tr>
<td class="org-right">17</td>
<td class="org-right">-0.413</td>
<td class="org-right">0.045</td>
<td class="org-right">0.118</td>
<td class="org-right">0.032</td>
</tr align="right">


<tr>
<td class="org-right">18</td>
<td class="org-right">-0.289</td>
<td class="org-right">0.042</td>
<td class="org-right">0.12</td>
<td class="org-right">0.031</td>
</tr align="right">


<tr>
<td class="org-right">19</td>
<td class="org-right">-0.412</td>
<td class="org-right">0.044</td>
<td class="org-right">0.087</td>
<td class="org-right">0.029</td>
</tr align="right">


<tr>
<td class="org-right">20</td>
<td class="org-right">-0.1</td>
<td class="org-right">0.038</td>
<td class="org-right">0.048</td>
<td class="org-right">0.028</td>
</tr align="right">


<tr>
<td class="org-right">21</td>
<td class="org-right">-0.942</td>
<td class="org-right">0.038</td>
<td class="org-right">0.575</td>
<td class="org-right">0.039</td>
</tr align="right">


<tr>
<td class="org-right">22</td>
<td class="org-right">-0.838</td>
<td class="org-right">0.041</td>
<td class="org-right">0.374</td>
<td class="org-right">0.038</td>
</tr align="right">


<tr>
<td class="org-right">23</td>
<td class="org-right">-0.8</td>
<td class="org-right">0.041</td>
<td class="org-right">0.349</td>
<td class="org-right">0.038</td>
</tr align="right">


<tr>
<td class="org-right">24</td>
<td class="org-right">-0.32</td>
<td class="org-right">0.043</td>
<td class="org-right">0.021</td>
<td class="org-right">0.026</td>
</tr align="right">


<tr>
<td class="org-right">25</td>
<td class="org-right">-0.768</td>
<td class="org-right">0.041</td>
<td class="org-right">0.331</td>
<td class="org-right">0.038</td>
</tr align="right">


<tr>
<td class="org-right">26</td>
<td class="org-right">-0.316</td>
<td class="org-right">0.044</td>
<td class="org-right">-0.017</td>
<td class="org-right">0.023</td>
</tr align="right">


<tr>
<td class="org-right">27</td>
<td class="org-right">-0.066</td>
<td class="org-right">0.033</td>
<td class="org-right">0.026</td>
<td class="org-right">0.025</td>
</tr align="right">


<tr>
<td class="org-right">28</td>
<td class="org-right">-0.267</td>
<td class="org-right">0.044</td>
<td class="org-right">0.095</td>
<td class="org-right">0.03</td>
</tr align="right">


<tr>
<td class="org-right">29</td>
<td class="org-right">-0.501</td>
<td class="org-right">0.044</td>
<td class="org-right">0.096</td>
<td class="org-right">0.031</td>
</tr align="right">


<tr>
<td class="org-right">30</td>
<td class="org-right">-0.046</td>
<td class="org-right">0.033</td>
<td class="org-right">0.137</td>
<td class="org-right">0.032</td>
</tr align="right">


<tr>
<td class="org-right">31</td>
<td class="org-right">-0.955</td>
<td class="org-right">0.038</td>
<td class="org-right">0.711</td>
<td class="org-right">0.041</td>
</tr align="right">


<tr>
<td class="org-right">32</td>
<td class="org-right">-0.933</td>
<td class="org-right">0.038</td>
<td class="org-right">0.643</td>
<td class="org-right">0.04</td>
</tr align="right">


<tr>
<td class="org-right">33</td>
<td class="org-right">-0.92</td>
<td class="org-right">0.039</td>
<td class="org-right">0.432</td>
<td class="org-right">0.039</td>
</tr align="right">


<tr>
<td class="org-right">34</td>
<td class="org-right">-0.425</td>
<td class="org-right">0.044</td>
<td class="org-right">-0.009</td>
<td class="org-right">0.024</td>
</tr align="right">


<tr>
<td class="org-right">35</td>
<td class="org-right">-0.51</td>
<td class="org-right">0.043</td>
<td class="org-right">0.065</td>
<td class="org-right">0.028</td>
</tr align="right">


<tr>
<td class="org-right">36</td>
<td class="org-right">-0.473</td>
<td class="org-right">0.043</td>
<td class="org-right">0.1</td>
<td class="org-right">0.03</td>
</tr align="right">


<tr>
<td class="org-right">37</td>
<td class="org-right">-0.135</td>
<td class="org-right">0.038</td>
<td class="org-right">0.008</td>
<td class="org-right">0.025</td>
</tr align="right">


<tr>
<td class="org-right">38</td>
<td class="org-right">-0.257</td>
<td class="org-right">0.04</td>
<td class="org-right">-0.005</td>
<td class="org-right">0.024</td>
</tr align="right">


<tr>
<td class="org-right">39</td>
<td class="org-right">-0.383</td>
<td class="org-right">0.044</td>
<td class="org-right">0.043</td>
<td class="org-right">0.027</td>
</tr align="right">


<tr>
<td class="org-right">40</td>
<td class="org-right">-0.125</td>
<td class="org-right">0.037</td>
<td class="org-right">0.167</td>
<td class="org-right">0.033</td>
</tr align="right">
</tbody>
</table>


## Log loss differences by response type

<table border="2" cellspacing="0" cellpadding="6" rules="groups" frame="hsides">


<colgroup>
<col  class="org-left" />

<col  class="org-right" />

<col  class="org-right" />

<col  class="org-right" />

<col  class="org-right" />

<col  class="org-right" />
</colgroup>
<thead>
<tr>
<th scope="col" class="org-left">res</th>
<th scope="col" class="org-right">item</th>
<th scope="col" class="org-right">pos\<sub>neg</sub></th>
<th scope="col" class="org-right">se\<sub>pn</sub></th>
<th scope="col" class="org-right">neg\<sub>pos</sub></th>
<th scope="col" class="org-right">se\<sub>np</sub></th>
</tr align="right">
</thead>

<tbody>
<tr>
<td class="org-left">incorrect</td>
<td class="org-right">1</td>
<td class="org-right">-0.083</td>
<td class="org-right">0.014</td>
<td class="org-right">-0.126</td>
<td class="org-right">0.018</td>
</tr align="right">


<tr>
<td class="org-left">correct</td>
<td class="org-right">1</td>
<td class="org-right">-0.970</td>
<td class="org-right">0.038</td>
<td class="org-right">0.726</td>
<td class="org-right">0.041</td>
</tr align="right">


<tr>
<td class="org-left">incorrect</td>
<td class="org-right">2</td>
<td class="org-right">-0.146</td>
<td class="org-right">0.019</td>
<td class="org-right">-0.318</td>
<td class="org-right">0.034</td>
</tr align="right">


<tr>
<td class="org-left">correct</td>
<td class="org-right">2</td>
<td class="org-right">-0.995</td>
<td class="org-right">0.037</td>
<td class="org-right">0.696</td>
<td class="org-right">0.037</td>
</tr align="right">


<tr>
<td class="org-left">incorrect</td>
<td class="org-right">3</td>
<td class="org-right">-0.164</td>
<td class="org-right">0.015</td>
<td class="org-right">-0.270</td>
<td class="org-right">0.026</td>
</tr align="right">


<tr>
<td class="org-left">correct</td>
<td class="org-right">3</td>
<td class="org-right">-1.051</td>
<td class="org-right">0.037</td>
<td class="org-right">0.620</td>
<td class="org-right">0.035</td>
</tr align="right">


<tr>
<td class="org-left">incorrect</td>
<td class="org-right">4</td>
<td class="org-right">-0.237</td>
<td class="org-right">0.009</td>
<td class="org-right">-0.398</td>
<td class="org-right">0.024</td>
</tr align="right">


<tr>
<td class="org-left">correct</td>
<td class="org-right">4</td>
<td class="org-right">-1.033</td>
<td class="org-right">0.037</td>
<td class="org-right">0.634</td>
<td class="org-right">0.035</td>
</tr align="right">


<tr>
<td class="org-left">incorrect</td>
<td class="org-right">5</td>
<td class="org-right">-0.178</td>
<td class="org-right">0.016</td>
<td class="org-right">-0.332</td>
<td class="org-right">0.025</td>
</tr align="right">


<tr>
<td class="org-left">correct</td>
<td class="org-right">5</td>
<td class="org-right">-1.012</td>
<td class="org-right">0.037</td>
<td class="org-right">0.650</td>
<td class="org-right">0.036</td>
</tr align="right">


<tr>
<td class="org-left">incorrect</td>
<td class="org-right">6</td>
<td class="org-right">-0.203</td>
<td class="org-right">0.010</td>
<td class="org-right">-0.318</td>
<td class="org-right">0.024</td>
</tr align="right">


<tr>
<td class="org-left">correct</td>
<td class="org-right">6</td>
<td class="org-right">-1.109</td>
<td class="org-right">0.035</td>
<td class="org-right">0.542</td>
<td class="org-right">0.031</td>
</tr align="right">


<tr>
<td class="org-left">incorrect</td>
<td class="org-right">7</td>
<td class="org-right">-0.200</td>
<td class="org-right">0.013</td>
<td class="org-right">-0.324</td>
<td class="org-right">0.027</td>
</tr align="right">


<tr>
<td class="org-left">correct</td>
<td class="org-right">7</td>
<td class="org-right">-1.109</td>
<td class="org-right">0.034</td>
<td class="org-right">0.541</td>
<td class="org-right">0.029</td>
</tr align="right">


<tr>
<td class="org-left">incorrect</td>
<td class="org-right">8</td>
<td class="org-right">-0.208</td>
<td class="org-right">0.012</td>
<td class="org-right">-0.292</td>
<td class="org-right">0.023</td>
</tr align="right">


<tr>
<td class="org-left">correct</td>
<td class="org-right">8</td>
<td class="org-right">-1.174</td>
<td class="org-right">0.034</td>
<td class="org-right">0.474</td>
<td class="org-right">0.029</td>
</tr align="right">


<tr>
<td class="org-left">incorrect</td>
<td class="org-right">9</td>
<td class="org-right">-0.107</td>
<td class="org-right">0.016</td>
<td class="org-right">-0.149</td>
<td class="org-right">0.023</td>
</tr align="right">


<tr>
<td class="org-left">correct</td>
<td class="org-right">9</td>
<td class="org-right">-1.317</td>
<td class="org-right">0.033</td>
<td class="org-right">0.348</td>
<td class="org-right">0.026</td>
</tr align="right">


<tr>
<td class="org-left">incorrect</td>
<td class="org-right">10</td>
<td class="org-right">-0.113</td>
<td class="org-right">0.016</td>
<td class="org-right">-0.161</td>
<td class="org-right">0.024</td>
</tr align="right">


<tr>
<td class="org-left">correct</td>
<td class="org-right">10</td>
<td class="org-right">-1.078</td>
<td class="org-right">0.038</td>
<td class="org-right">0.589</td>
<td class="org-right">0.036</td>
</tr align="right">


<tr>
<td class="org-left">incorrect</td>
<td class="org-right">11</td>
<td class="org-right">-0.167</td>
<td class="org-right">0.014</td>
<td class="org-right">-0.221</td>
<td class="org-right">0.022</td>
</tr align="right">


<tr>
<td class="org-left">correct</td>
<td class="org-right">11</td>
<td class="org-right">-1.278</td>
<td class="org-right">0.033</td>
<td class="org-right">0.372</td>
<td class="org-right">0.027</td>
</tr align="right">


<tr>
<td class="org-left">incorrect</td>
<td class="org-right">12</td>
<td class="org-right">-0.106</td>
<td class="org-right">0.018</td>
<td class="org-right">-0.154</td>
<td class="org-right">0.025</td>
</tr align="right">


<tr>
<td class="org-left">correct</td>
<td class="org-right">12</td>
<td class="org-right">-1.325</td>
<td class="org-right">0.032</td>
<td class="org-right">0.305</td>
<td class="org-right">0.018</td>
</tr align="right">


<tr>
<td class="org-left">incorrect</td>
<td class="org-right">13</td>
<td class="org-right">-0.100</td>
<td class="org-right">0.018</td>
<td class="org-right">-0.155</td>
<td class="org-right">0.026</td>
</tr align="right">


<tr>
<td class="org-left">correct</td>
<td class="org-right">13</td>
<td class="org-right">-1.170</td>
<td class="org-right">0.034</td>
<td class="org-right">0.458</td>
<td class="org-right">0.026</td>
</tr align="right">


<tr>
<td class="org-left">incorrect</td>
<td class="org-right">14</td>
<td class="org-right">-0.097</td>
<td class="org-right">0.016</td>
<td class="org-right">-0.143</td>
<td class="org-right">0.024</td>
</tr align="right">


<tr>
<td class="org-left">correct</td>
<td class="org-right">14</td>
<td class="org-right">-1.329</td>
<td class="org-right">0.033</td>
<td class="org-right">0.335</td>
<td class="org-right">0.022</td>
</tr align="right">


<tr>
<td class="org-left">incorrect</td>
<td class="org-right">15</td>
<td class="org-right">-0.179</td>
<td class="org-right">0.014</td>
<td class="org-right">-0.253</td>
<td class="org-right">0.024</td>
</tr align="right">


<tr>
<td class="org-left">correct</td>
<td class="org-right">15</td>
<td class="org-right">-1.208</td>
<td class="org-right">0.034</td>
<td class="org-right">0.444</td>
<td class="org-right">0.025</td>
</tr align="right">


<tr>
<td class="org-left">incorrect</td>
<td class="org-right">16</td>
<td class="org-right">-0.030</td>
<td class="org-right">0.019</td>
<td class="org-right">-0.065</td>
<td class="org-right">0.024</td>
</tr align="right">


<tr>
<td class="org-left">correct</td>
<td class="org-right">16</td>
<td class="org-right">-1.466</td>
<td class="org-right">0.032</td>
<td class="org-right">0.184</td>
<td class="org-right">0.016</td>
</tr align="right">


<tr>
<td class="org-left">incorrect</td>
<td class="org-right">17</td>
<td class="org-right">-0.018</td>
<td class="org-right">0.020</td>
<td class="org-right">-0.055</td>
<td class="org-right">0.025</td>
</tr align="right">


<tr>
<td class="org-left">correct</td>
<td class="org-right">17</td>
<td class="org-right">-1.198</td>
<td class="org-right">0.042</td>
<td class="org-right">0.461</td>
<td class="org-right">0.032</td>
</tr align="right">


<tr>
<td class="org-left">incorrect</td>
<td class="org-right">18</td>
<td class="org-right">0.015</td>
<td class="org-right">0.021</td>
<td class="org-right">-0.019</td>
<td class="org-right">0.025</td>
</tr align="right">


<tr>
<td class="org-left">correct</td>
<td class="org-right">18</td>
<td class="org-right">-1.106</td>
<td class="org-right">0.040</td>
<td class="org-right">0.495</td>
<td class="org-right">0.034</td>
</tr align="right">


<tr>
<td class="org-left">incorrect</td>
<td class="org-right">19</td>
<td class="org-right">-0.012</td>
<td class="org-right">0.023</td>
<td class="org-right">-0.067</td>
<td class="org-right">0.027</td>
</tr align="right">


<tr>
<td class="org-left">correct</td>
<td class="org-right">19</td>
<td class="org-right">-1.223</td>
<td class="org-right">0.030</td>
<td class="org-right">0.401</td>
<td class="org-right">0.023</td>
</tr align="right">


<tr>
<td class="org-left">incorrect</td>
<td class="org-right">20</td>
<td class="org-right">0.036</td>
<td class="org-right">0.021</td>
<td class="org-right">0.013</td>
<td class="org-right">0.024</td>
</tr align="right">


<tr>
<td class="org-left">correct</td>
<td class="org-right">20</td>
<td class="org-right">-1.444</td>
<td class="org-right">0.058</td>
<td class="org-right">0.395</td>
<td class="org-right">0.051</td>
</tr align="right">


<tr>
<td class="org-left">incorrect</td>
<td class="org-right">21</td>
<td class="org-right">-0.193</td>
<td class="org-right">0.011</td>
<td class="org-right">-0.401</td>
<td class="org-right">0.032</td>
</tr align="right">


<tr>
<td class="org-left">correct</td>
<td class="org-right">21</td>
<td class="org-right">-1.010</td>
<td class="org-right">0.036</td>
<td class="org-right">0.664</td>
<td class="org-right">0.035</td>
</tr align="right">


<tr>
<td class="org-left">incorrect</td>
<td class="org-right">22</td>
<td class="org-right">-0.132</td>
<td class="org-right">0.016</td>
<td class="org-right">-0.231</td>
<td class="org-right">0.030</td>
</tr align="right">


<tr>
<td class="org-left">correct</td>
<td class="org-right">22</td>
<td class="org-right">-1.079</td>
<td class="org-right">0.035</td>
<td class="org-right">0.581</td>
<td class="org-right">0.031</td>
</tr align="right">


<tr>
<td class="org-left">incorrect</td>
<td class="org-right">23</td>
<td class="org-right">-0.114</td>
<td class="org-right">0.019</td>
<td class="org-right">-0.205</td>
<td class="org-right">0.029</td>
</tr align="right">


<tr>
<td class="org-left">correct</td>
<td class="org-right">23</td>
<td class="org-right">-1.081</td>
<td class="org-right">0.033</td>
<td class="org-right">0.576</td>
<td class="org-right">0.031</td>
</tr align="right">


<tr>
<td class="org-left">incorrect</td>
<td class="org-right">24</td>
<td class="org-right">-0.020</td>
<td class="org-right">0.019</td>
<td class="org-right">-0.055</td>
<td class="org-right">0.025</td>
</tr align="right">


<tr>
<td class="org-left">correct</td>
<td class="org-right">24</td>
<td class="org-right">-1.341</td>
<td class="org-right">0.036</td>
<td class="org-right">0.278</td>
<td class="org-right">0.023</td>
</tr align="right">


<tr>
<td class="org-left">incorrect</td>
<td class="org-right">25</td>
<td class="org-right">-0.105</td>
<td class="org-right">0.019</td>
<td class="org-right">-0.185</td>
<td class="org-right">0.026</td>
</tr align="right">


<tr>
<td class="org-left">correct</td>
<td class="org-right">25</td>
<td class="org-right">-1.079</td>
<td class="org-right">0.033</td>
<td class="org-right">0.573</td>
<td class="org-right">0.033</td>
</tr align="right">


<tr>
<td class="org-left">incorrect</td>
<td class="org-right">26</td>
<td class="org-right">-0.029</td>
<td class="org-right">0.018</td>
<td class="org-right">-0.060</td>
<td class="org-right">0.024</td>
</tr align="right">


<tr>
<td class="org-left">correct</td>
<td class="org-right">26</td>
<td class="org-right">-1.499</td>
<td class="org-right">0.037</td>
<td class="org-right">0.161</td>
<td class="org-right">0.016</td>
</tr align="right">


<tr>
<td class="org-left">incorrect</td>
<td class="org-right">27</td>
<td class="org-right">0.043</td>
<td class="org-right">0.022</td>
<td class="org-right">0.009</td>
<td class="org-right">0.026</td>
</tr align="right">


<tr>
<td class="org-left">correct</td>
<td class="org-right">27</td>
<td class="org-right">-1.326</td>
<td class="org-right">0.029</td>
<td class="org-right">0.217</td>
<td class="org-right">0.021</td>
</tr align="right">


<tr>
<td class="org-left">incorrect</td>
<td class="org-right">28</td>
<td class="org-right">0.026</td>
<td class="org-right">0.023</td>
<td class="org-right">-0.013</td>
<td class="org-right">0.026</td>
</tr align="right">


<tr>
<td class="org-left">correct</td>
<td class="org-right">28</td>
<td class="org-right">-1.242</td>
<td class="org-right">0.040</td>
<td class="org-right">0.456</td>
<td class="org-right">0.033</td>
</tr align="right">


<tr>
<td class="org-left">incorrect</td>
<td class="org-right">29</td>
<td class="org-right">-0.055</td>
<td class="org-right">0.019</td>
<td class="org-right">-0.096</td>
<td class="org-right">0.025</td>
</tr align="right">


<tr>
<td class="org-left">correct</td>
<td class="org-right">29</td>
<td class="org-right">-1.233</td>
<td class="org-right">0.031</td>
<td class="org-right">0.412</td>
<td class="org-right">0.030</td>
</tr align="right">


<tr>
<td class="org-left">incorrect</td>
<td class="org-right">30</td>
<td class="org-right">0.076</td>
<td class="org-right">0.023</td>
<td class="org-right">0.045</td>
<td class="org-right">0.025</td>
</tr align="right">


<tr>
<td class="org-left">correct</td>
<td class="org-right">30</td>
<td class="org-right">-0.907</td>
<td class="org-right">0.044</td>
<td class="org-right">0.786</td>
<td class="org-right">0.045</td>
</tr align="right">


<tr>
<td class="org-left">incorrect</td>
<td class="org-right">31</td>
<td class="org-right">-0.141</td>
<td class="org-right">0.006</td>
<td class="org-right">-0.259</td>
<td class="org-right">0.017</td>
</tr align="right">


<tr>
<td class="org-left">correct</td>
<td class="org-right">31</td>
<td class="org-right">-0.969</td>
<td class="org-right">0.038</td>
<td class="org-right">0.727</td>
<td class="org-right">0.040</td>
</tr align="right">


<tr>
<td class="org-left">incorrect</td>
<td class="org-right">32</td>
<td class="org-right">-0.064</td>
<td class="org-right">0.021</td>
<td class="org-right">-0.249</td>
<td class="org-right">0.038</td>
</tr align="right">


<tr>
<td class="org-left">correct</td>
<td class="org-right">32</td>
<td class="org-right">-0.989</td>
<td class="org-right">0.036</td>
<td class="org-right">0.700</td>
<td class="org-right">0.037</td>
</tr align="right">


<tr>
<td class="org-left">incorrect</td>
<td class="org-right">33</td>
<td class="org-right">-0.218</td>
<td class="org-right">0.013</td>
<td class="org-right">-0.352</td>
<td class="org-right">0.025</td>
</tr align="right">


<tr>
<td class="org-left">correct</td>
<td class="org-right">33</td>
<td class="org-right">-1.065</td>
<td class="org-right">0.036</td>
<td class="org-right">0.594</td>
<td class="org-right">0.033</td>
</tr align="right">


<tr>
<td class="org-left">incorrect</td>
<td class="org-right">34</td>
<td class="org-right">-0.057</td>
<td class="org-right">0.018</td>
<td class="org-right">-0.094</td>
<td class="org-right">0.024</td>
</tr align="right">


<tr>
<td class="org-left">correct</td>
<td class="org-right">34</td>
<td class="org-right">-1.414</td>
<td class="org-right">0.032</td>
<td class="org-right">0.221</td>
<td class="org-right">0.016</td>
</tr align="right">


<tr>
<td class="org-left">incorrect</td>
<td class="org-right">35</td>
<td class="org-right">-0.062</td>
<td class="org-right">0.019</td>
<td class="org-right">-0.113</td>
<td class="org-right">0.025</td>
</tr align="right">


<tr>
<td class="org-left">correct</td>
<td class="org-right">35</td>
<td class="org-right">-1.247</td>
<td class="org-right">0.029</td>
<td class="org-right">0.359</td>
<td class="org-right">0.022</td>
</tr align="right">


<tr>
<td class="org-left">incorrect</td>
<td class="org-right">36</td>
<td class="org-right">-0.037</td>
<td class="org-right">0.020</td>
<td class="org-right">-0.089</td>
<td class="org-right">0.026</td>
</tr align="right">


<tr>
<td class="org-left">correct</td>
<td class="org-right">36</td>
<td class="org-right">-1.201</td>
<td class="org-right">0.030</td>
<td class="org-right">0.415</td>
<td class="org-right">0.025</td>
</tr align="right">


<tr>
<td class="org-left">incorrect</td>
<td class="org-right">37</td>
<td class="org-right">0.022</td>
<td class="org-right">0.020</td>
<td class="org-right">-0.007</td>
<td class="org-right">0.025</td>
</tr align="right">


<tr>
<td class="org-left">correct</td>
<td class="org-right">37</td>
<td class="org-right">-1.552</td>
<td class="org-right">0.045</td>
<td class="org-right">0.146</td>
<td class="org-right">0.021</td>
</tr align="right">


<tr>
<td class="org-left">incorrect</td>
<td class="org-right">38</td>
<td class="org-right">-0.012</td>
<td class="org-right">0.020</td>
<td class="org-right">-0.045</td>
<td class="org-right">0.025</td>
</tr align="right">


<tr>
<td class="org-left">correct</td>
<td class="org-right">38</td>
<td class="org-right">-1.413</td>
<td class="org-right">0.024</td>
<td class="org-right">0.183</td>
<td class="org-right">0.014</td>
</tr align="right">


<tr>
<td class="org-left">incorrect</td>
<td class="org-right">39</td>
<td class="org-right">-0.028</td>
<td class="org-right">0.020</td>
<td class="org-right">-0.066</td>
<td class="org-right">0.025</td>
</tr align="right">


<tr>
<td class="org-left">correct</td>
<td class="org-right">39</td>
<td class="org-right">-1.318</td>
<td class="org-right">0.035</td>
<td class="org-right">0.331</td>
<td class="org-right">0.024</td>
</tr align="right">


<tr>
<td class="org-left">incorrect</td>
<td class="org-right">40</td>
<td class="org-right">0.072</td>
<td class="org-right">0.023</td>
<td class="org-right">0.038</td>
<td class="org-right">0.026</td>
</tr align="right">


<tr>
<td class="org-left">correct</td>
<td class="org-right">40</td>
<td class="org-right">-0.961</td>
<td class="org-right">0.041</td>
<td class="org-right">0.712</td>
<td class="org-right">0.039</td>
</tr align="right">
</tbody>
</table>

