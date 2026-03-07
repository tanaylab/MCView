# Create server logic for copy button notifications

Create server logic for copy button notifications

## Usage

``` r
clipboard_copy_button_server(
  input,
  id,
  data_reactive,
  globals = NULL,
  message_template = "Copied {count} items to clipboard"
)
```

## Arguments

- input:

  Shiny input

- id:

  Button ID

- data_reactive:

  Reactive expression that returns the data

- globals:

  Global reactive values (optional, for internal clipboard)

- message_template:

  Template for success message with count placeholder

## Value

Observer for button clicks
