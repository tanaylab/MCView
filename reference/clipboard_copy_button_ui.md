# Create a reactive copy button using rclipboard

Create a reactive copy button using rclipboard

## Usage

``` r
clipboard_copy_button_ui(
  ns,
  id,
  data_reactive,
  label = "Copy to Clipboard",
  style = "background-color: #17a2b8; color: white; border: none;",
  tooltip = "Copy to system clipboard",
  disabled_label = NULL
)
```

## Arguments

- ns:

  Namespace function

- id:

  Button ID

- data_reactive:

  Reactive expression that returns the data to copy

- label:

  Button label

- style:

  Button styling

- tooltip:

  Tooltip text

- disabled_label:

  Label when button is disabled

## Value

UI output for the copy button
