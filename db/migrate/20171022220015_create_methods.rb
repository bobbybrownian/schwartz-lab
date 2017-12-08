class CreateMethods < ActiveRecord::Migration[5.1]
  def change
    create_table :labmethods do |t|

        t.string :title, null: false
        t.string :pdf_url, null: false
        t.belongs_to :methodtype, null: false

        t.timestamps null: false
    end
  end
end
