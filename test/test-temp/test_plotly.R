

p <- plot_ly(data = iris, x = ~Sepal.Length, y = ~Petal.Length,
             marker = list(size = 10,
                           color = 'rgba(255, 182, 193, .9)',
                           line = list(color = 'rgba(152, 0, 0, .8)',
                                       width = 2))) %>%
    layout(title = 'Styled Scatter',
           yaxis = list(zeroline = FALSE),
           xaxis = list(zeroline = FALSE))

# Create a shareable link to your chart
# Set up API credentials: https://plot.ly/r/getting-started
# chart_link = api_create(p, filename="scatter-styled")
# chart_link
p

{ggplot(dnew, aes(t, y)) + geom_line() + geom_point(aes(color = wf, shape = wf))} %>% ggplotly
